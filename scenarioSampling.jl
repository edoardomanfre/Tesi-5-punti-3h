#using Distributions
#using StatsPlots
#using CSV
#using Printf
#using Clustering
#using Plots
#plotly()



# helping functions
#--------------------------------------------------------------
# normalize data
function normalize(data)
    # data mxn matrix, m=stages, n=samples
    dataMean = mean(data, dims = 2)
    dataStd = std(data, dims = 2)

    newData = zeros(size(data))

    for i = 1:size(data)[1]
        newData[i, :] = (data[i, :] .- dataMean[i]) ./ dataStd[i]
    end
    return NormData(newData, dataMean, dataStd)
end

# find correlation
function correlation(v1, v2)                                        # v1= inflow, v2 = prezzo
    correlation = zeros(size(v1)[2])                                # Creo vettore correlazione di dimensione 30
    for s = 1:size(v1)[2]                                           #s=1:30 anni
        correlation[s] = cor(v1[:, s], v2[:, s])                    # Per ogni anno, si trova la correlazione settimanala tre inflow e prezzo
    end
    return mean(correlation)                                        # Trova la corrleazione media tra i 30 anni
end

# sample random value from distribution
function sample_from_distribution(data, distribution, nSamples)
    if distribution == "Normal"
        dist = MvNormal(vec(mean(data, dims = 1)), cov(data, dims = 1))
    elseif distribution == "LogNormal"
        dist = MvLogNormal(vec(mean(data, dims = 1)), cov(data, dims = 1))
    else
        error("Did not match distribution: use \"Normal\" or \"LogNormal\"")
    end
    #Random.seed!(485)
    samples = rand(dist, nSamples)
    return samples
end

# replace all negative values 
function replace_negatives(
    data_norm,
    distribution,
    trajectories,
    trajectories_corrected,
    dataI,
    dataP,
    corr_I,
    corr_IP,
    seriesim
)
    negatives = findall(x -> x < 0, trajectories_corrected)
    # Draw new samples for negative values and adjust for correlation, mean and std
    while length(negatives) > 0
        for n in negatives
            s = sample_from_distribution(data_norm, distribution, 1)
            if n[1] > 1
                trajectories[n[1], 1, n[3]] = s[1] + trajectories[n[1]-1, 1, n[3]] * corr_I
                trajectories[n[1], 2, n[3]] =
                    s[2] + (trajectories[n[1], 1, n[3]] - s[1]) * corr_IP
            elseif n[1] ==1 && n[3] !=1 
                if seriesim
                    trajectories[n[1], 1, n[3]] =
                    s[1] .+ trajectories[end, 1, n[3]-1] .* corr_I
                    trajectories[n[1], 2, n[3]] =
                        s[2] + (trajectories[n[1], 1, n[3]] - s[1]) * corr_IP
                else
                    trajectories[n[1], :, n[3]] = s
                end
            else
                trajectories[n[1], :, n[3]] = s
            end
            trajectories_corrected[n[1], 1, n[3]] =
                (trajectories[n[1], 1, n[3]] * dataI.std[n[1]]) .+ dataI.mean[n[1]]
            trajectories_corrected[n[1], 2, n[3]] =
                (trajectories[n[1], 2, n[3]] * dataP.std[n[1]]) .+ dataP.mean[n[1]]
        end
        negatives = findall(x -> x < 0, trajectories_corrected)
    end


    # Adjust to compensate for change in mean
    for t = 1:size(dataI.data)[1]
        trajectories_corrected[t, 1, :] =
            (trajectories_corrected[t, 1, :] ./ mean(trajectories_corrected[t, 1, :])) .*
            dataI.mean[t]
        trajectories_corrected[t, 2, :] =
            (trajectories_corrected[t, 2, :] ./ mean(trajectories_corrected[t, 2, :])) .*
            dataP.mean[t]
    end

    return trajectories_corrected, trajectories
end

#sort states in markov model by inflow 
function sort_clusters_by_states(KmeansClusters)
    SortedKmeansClusters = Vector{SortedClusters}(undef, length(KmeansClusters))
    for t = 1:length(KmeansClusters)
        centers_new = zeros(size(KmeansClusters[t].centers))
        counts_new = zeros(size(KmeansClusters[t].centers)[2])
        assignments_new = zeros(length(KmeansClusters[t].assignments))
        for state = 1:size(KmeansClusters[t].centers)[2]
            temp = findmax(KmeansClusters[t].centers[1, :])
            loc = temp[2]
            centers_new[1, state] = temp[1]
            centers_new[2, state] = KmeansClusters[t].centers[2, loc]
            counts_new[state] = KmeansClusters[t].counts[loc]
            KmeansClusters[t].centers[:, loc] .= -10000
            for scen = 1:length(KmeansClusters[t].assignments)
                if KmeansClusters[t].assignments[scen] == loc
                    assignments_new[scen] = state
                end
            end
        end
        SortedKmeansClusters[t] = SortedClusters(centers_new, counts_new, assignments_new)
    end
    return SortedKmeansClusters
end

# Generate stochastic model for use in SDP
#-----------------------------------------------------------------

# make Markov model
function samplingAlg(inflow, price, InputParameters, seednum = 132, serial=true)
    #Sample trajectories
    @unpack (
        NStage,
        NStates,
        NSamples,
    ) = InputParameters

    trajectories, trajectories_normalized, dataI, dataP =
        make_trajectories(inflow, price,  NSamples, "Normal", seednum, serial)

    #R = [kmeans(trajectories[t,:,:], NStates;maxiter=200) for t =1:T]
    KmeansClusters =
        [kmeans(trajectories_normalized[t, :, :], NStates; maxiter = 200) for t = 1:NStage]

    KmeansClusters = sort_clusters_by_states(KmeansClusters)

    @assert(
        length(KmeansClusters[1].assignments) == size(trajectories)[3],
        "All trajectories are not assigned!"
    )

    transitions = count_transitions(KmeansClusters)
    startProbability = KmeansClusters[1].counts ./ sum(KmeansClusters[1].counts)
    transitionProb = calculate_probabilities(KmeansClusters, transitions)
    states = calculate_states_from_normalized_data(KmeansClusters, dataI, dataP)
    envProbability = zeros(NStage, NStates, 2, 2)
    envProbability[:, :, 1, 1] .= 1
    ScenarioLattice = lattice(states, transitionProb, envProbability)

    return samplingData(
        trajectories,
        trajectories_normalized,
        dataI,
        dataP,
        KmeansClusters,
        transitions,
        startProbability,
        transitionProb,
        states,
        ScenarioLattice,
    )
end

# make autocorrelate scenario trajectories
function make_trajectories(inflow, price, NSamples = 10000, distribution = "Normal",  seednum=123, seriesim=true)
    # calculte correlation
    corr_IP = correlation(inflow, price) #average over all scen of the correlation between the inflow and price
    corr_I = correlation(inflow, vcat(inflow[2:end, :], inflow[1, :]')) #average over all scen of the correlation between the inflow in t and t+1 for each year

    #Normalize data and calculate parameters
    dataI = normalize(inflow)
    dataP = normalize(price)
    data_norm =
        [reshape(dataI.data, (length(dataI.data))) reshape(dataP.data, (length(dataP.data)))]

    NVariables = size(data_norm)[2]
    trajectories = zeros(size(dataI.data)[1], NVariables, NSamples)
    trajectories_corrected = zeros(size(dataI.data)[1], NVariables, NSamples)
    Random.seed!(seednum) 
    for t = 1:size(dataI.data)[1]
        trajectories[t, :, :] = sample_from_distribution(data_norm, distribution, NSamples)
        if t > 1
            temp = trajectories[t, 1, :]
            trajectories[t, 1, :] =
                trajectories[t, 1, :] .+ trajectories[t-1, 1, :] .* corr_I                              # Inflo
            trajectories[t, 2, :] =
                trajectories[t, 2, :] .+ (trajectories[t, 1, :] .- temp) .* corr_IP 
        end

    end

    if seriesim
        temp = trajectories[1, 1,2:end]                 # Inflow della settimana 1 per tutte le 9999 trajectories (non si considera la prima)
        trajectories[1, 1, 2:end] =
            trajectories[1, 1, 2:end] .+ trajectories[end, 1, 1:end-1] .* corr_I
        trajectories[1, 2, 2:end] =
            trajectories[1, 2, 2:end] .+ ( trajectories[1, 1, 2:end] .- temp) .* corr_IP 
    end
    
    for t = 1:size(dataI.data)[1]
        trajectories_corrected[t, 1, :] =
        (trajectories[t, 1, :] .* dataI.std[t]) .+ dataI.mean[t]
        trajectories_corrected[t, 2, :] =
        (trajectories[t, 2, :] .* dataP.std[t]) .+ dataP.mean[t]
    end

    trajectories_corrected, trajectories = replace_negatives(               # Rimpiazza i valori negativi con valori positivi o nulli
        data_norm,
        distribution,
        trajectories,
        trajectories_corrected,
        dataI,
        dataP,
        corr_I,
        corr_IP,
        seriesim,
    )


    #print some key parameters desctibing the input data and the sampled data
    println(" ")
    @printf(
        "\nInflow Input: \t\tmean: %.2f \tstd: %.2f \tmax: %.2f \tmin: %.2f\n",
        mean(inflow),
        std(inflow),
        findmax(inflow)[1],
        findmin(inflow)[1]
    )
    @printf(
        "Inflow Samples: \tmean: %.2f \tstd: %.2f \tmax: %.2f \tmin: %.2f\n",
        mean(trajectories_corrected[:, 1, :]),
        std(trajectories_corrected[:, 1, :])[1],
        findmax(trajectories_corrected[:, 1, :])[1],
        findmin(trajectories_corrected[:, 1, :])[1]
    )
    @printf(
        "\nPrice Input: \tmean: %.2f \tstd: %.2f \tmax: %.2f \tmin: %.2f\n",
        mean(price),
        std(price),
        findmax(price)[1],
        findmin(price)[1]
    )
    @printf(
        "Price Samples: \tmean: %.2f \tstd: %.2f \tmax: %.2f \tmin: %.2f\n",
        mean(trajectories_corrected[:, 2, :]),
        std(trajectories_corrected[:, 2, :]),
        findmax(trajectories_corrected[:, 2, :])[1],
        findmin(trajectories_corrected[:, 2, :])[1]
    )

    corr_IP = correlation(inflow, price) #average over all scen of the correlation between the inflow and price
    corr_I = correlation(inflow, vcat(inflow[2:end, :], inflow[1, :]'))
    corr_IP_sim = correlation(trajectories[:, 1, :], trajectories[:, 2, :]) #average over all scen of the correlation between the inflow and price
    corr_I_sim = correlation(
        trajectories[:, 1, :],
        vcat(trajectories[2:end, 1, :], trajectories[1, 1, :]'),
    )

    @printf(
        "\ncorrelation inflow-Price: \tinput: \t%.2f samples: \t%.2f\n",
        corr_IP,
        corr_IP_sim
    )
    @printf("correlation inflow: \t\tinput: \t%.2f samples: \t%.2f\n", corr_I, corr_I_sim)

    return trajectories_corrected, trajectories, dataI, dataP
end

# calculate transition probabilities
function count_transitions(R_norm)
    # transitions t gives transitions from t-1 to t (from T->1, 1->2,...,T-1->T)
    NTransitions = zeros(length(R_norm), length(R_norm[1].counts), length(R_norm[1].counts))                                # Creo vettore 52x5x5
    for t = 1:length(R_norm)                                                                                                # Inizia un ciclo per ogni settimana t=1:52
        if t == 1                                                                                                           # Se sono alla prima settimana, nFrom= R_norm[52].
            nFrom = R_norm[end].counts                                                                                      # Cluster della settimana 52 [1128,2388,2057,2219,2208]
            nTo = R_norm[end].counts
        else
            nFrom = R_norm[t-1].counts                                                                                      # Se siamo alle settimane successive, KmeansClusters[t-1].counts -> vettore di 5 valori
            nTo = R_norm[t-1].counts
        end

        for n = 1:length(R_norm[t].assignments)                                                                             # Inizia ciclo da n =1:10000 (per tutte le combinazioni)
            counted = false
            for to = 1:length(nTo)                                                                                          # per to=1:5
                for from = 1:length(nFrom)                                                                                  # Per from=1:5
                    if t == 1                                                                                               # Se sono alla prima settimana t=1
                        if R_norm[end].assignments[n] == from &&                                                            
                           R_norm[t].assignments[n] == to
                            NTransitions[t, from, to] += 1                                                                  # Conto quante combinazioni transitano da uno Stato N della settimana 52 allo stato T della settimana 1
                            counted = true
                            if counted
                                break
                            end
                        end
                    else
                        if R_norm[t-1].assignments[n] == from &&
                           R_norm[t].assignments[n] == to
                            NTransitions[t, from, to] += 1                                                                  # Conto quante combinazioni transitano da uno Stato N della settimana t allo stato T della settimana t+1
                            counted = true
                            if counted
                                break
                            end
                        end
                    end
                    if counted
                        break
                    end
                end
            end
        end
    end

    @assert(
        sum(NTransitions[1, :, :]) == length(R_norm[1].assignments),
        "Not all trajectories are accounted for."
    )
    @assert(
        sum(NTransitions[1, :, 1]) == sum(NTransitions[1+1, 1, :]),
        "In and out of first state does not add up."
    )
    @assert(
        sum(NTransitions[end-1, :, end]) == sum(NTransitions[end, end, :]),
        "In and out of last state does not add up."
    )

    return NTransitions
end

function calculate_probabilities(R_norm, transitions)
    #probabilities = zeros(size(transitions))
    probabilities = []
    for t = 1:size(transitions)[1]
        if t == 1
            append!(probabilities, [transitions[t, :, :] ./ R_norm[end].counts])
            #probabilities[t,:,:] .= transitions[t,:,:] ./ R_norm[end].counts
        else
            #probabilities[t,:,:] .= transitions[t,:,:] ./ R_norm[t-1].counts
            append!(probabilities, [transitions[t, :, :] ./ R_norm[t-1].counts])
        end
    end

   # t = rand(2:size(probabilities)[1])
   # nfrom = rand(1:size(probabilities[1])[1])
   # nto = rand(1:size(probabilities[1])[2])
    #for t = 2:size(probabilities)[1]
    #@assert(all(sum(probabilities[t], dims=2)) == 1.0, "\nSum of prob. not 1, but")
    #println("\n Sum probabilities: ", sum(probabilities[t],dims=2))
    #end

    return probabilities
end

# adjust back from normalized data
function calculate_states_from_normalized_data(KmeansClusters, dataI, dataP)
    states = []
    for t = 1:size(KmeansClusters)[1]
        stateVariables =
            [KmeansClusters[t].centers[1, :] .* dataI.std[t] .+ dataI.mean[t] KmeansClusters[t].centers[
                2,
                :,
            ] .* dataP.std[t] .+ dataP.mean[t]]
        append!(states, [stateVariables])
    end
    return states
end


# Make scenarios for end-simulation
#-----------------------------------------------------------------------

# draw from original trajecories    sorts 100 scenarios among all
function drawScenForSim(NSamples, NSimScen, trajectories, KmeansClusters)
    println("Drawing new scenarios..")
    Random.seed!(246)
    scen = rand(1:NSamples, NSimScen)
    scenarios = []
    scenStates = []
    for i in scen
        append!(scenarios, [trajectories[:, :, i]])
        append!(scenStates, [[Int(KmeansClusters[t].assignments[i]) for t = 1:NStage]])
    end
    return SimScenData(scenarios, scenStates)
end

# draw scenarios from new trajectories
function drawOutOfSampleScenForSim(priceYear, NSimScen, KmeansClusters, seednum=5236)
    println(string("Drawing new out-of-sample scenarios based on..", priceYear))
    
    Random.seed!(seednum)
    trajectories, trajectories_normalized, dataI, dataP =
        make_trajectories(inflow, price, NSimScen, "Normal", seednum, true)

    if length(priceYear)> 0
        println(string("\nScaling sampled prices to year ",priceYear ))
        priceFile = string("meanPrice_", priceYear, ".csv")
        meanPrice = Matrix{Float64}(
            CSV.read(joinpath(DataPath, priceFile), DataFrame, header = false),
        ) #EUR/MWh
        
        #newPricedata =  NormData(scenLData.dataP.data, meanPrice, scenLData.dataP.std);
        newTrajectories = zeros(size(trajectories));

        for s = 1:NSimScen
            for t= 1:size(trajectories)[1]
                newTrajectories[t,1,s] = trajectories[t,1,s]
                newTrajectories[t,2,s] = trajectories[t,2,s] .* (meanPrice[t] ./ dataP.mean[t])
            end
        end
    else
        newTrajectories= trajectories
    end


    scenarios = []
    scenStates = []
    for scen = 1:size(trajectories)[3]
        #println("Scen: ", scen)
        selectedStates = []
        min_diff_i = []
        scen_temp = zeros(size(trajectories)[1], 2)
        for t = 1:size(trajectories)[1]
            temp_i = abs.(newTrajectories[t, 1, scen] .- KmeansClusters[t][:, 1]) ./ mean(newTrajectories[t, 1,:])
            temp_p = abs.(newTrajectories[t, 2, scen] .- KmeansClusters[t][:, 2])  ./ mean(newTrajectories[t, 2,:])
            #temp_i = abs.(trajectories_normalized[t, 1, scen] .- KmeansClusters[t].centers[1, :])
            #temp_p = abs.(trajectories_normalized[t, 2, scen] .- KmeansClusters[t].centers[2, :])
            #temp_i = abs.(trajectories[t, 1, scen] .- states[t][:, 1]) ./ trajectories[t, 1, scen]
            #temp_p = abs.(trajectories[t, 2, scen] .- states[t][:, 2]) ./ trajectories[t, 2, scen]
            #append!(min_diff_i, findmin(temp_i)[2][1])
            append!(selectedStates, findmin(sqrt.(temp_i.^2 .+ temp_p.^2))[2][1])
            scen_temp[t, 1] = newTrajectories[t, 1, scen]
            scen_temp[t, 2] = newTrajectories[t, 2, scen]
        end
        append!(scenarios, [scen_temp])
        append!(scenStates, [selectedStates])
    end
    return SimScenData(scenarios, scenStates)
end

# use historical input data as scenarios 
function scenFromInput(inflow, price, states)
    scenarios = []
    scenStates = []
    for scen = 1:size(inflow)[2]
        #println("Scen: ", scen)
        selectedStates = []
        min_diff_i = []
        scen_temp = zeros(size(inflow)[1], 2)
        for t = 1:size(inflow)[1]
            temp_i = abs.(inflow[t, scen] .- states[t][:, 1]) ./ inflow[t, scen]
            temp_p = abs.(price[t, scen] .- states[t][:, 2]) ./ price[t, scen]
            append!(min_diff_i, findmin(temp_i)[2][1])
            append!(selectedStates, findmin(temp_i .+ temp_p)[2][1])
            scen_temp[t, 1] = inflow[t, scen]
            scen_temp[t, 2] = price[t, scen]
            #if findmin(temp_i)[2][1] != findmin(temp_i .+ temp_p)[2][1]
            #    @printf("%0.2f : %0.2f \n", findmin(temp_i)[1][1], temp_i[findmin(temp_i .+ temp_p)[2][1]])
            #end

        end
        append!(scenarios, [scen_temp])
        append!(scenStates, [selectedStates])
    end
    NSimScen = size(inflow)[2]

    return SimScenData(scenarios, scenStates), NSimScen
end

# adjust stochastic input
#--------------------------------------------------------------------------

# extend scenario lattice (markov model) to include environmetnal state
function extend_scenario_lattice(scenLData, envData)
    tracker = Calculate_envConst_prob(envData, scenLData, HY)

    scenLData.ScenarioLattice.envProbability[
        envData.firstAct+1:envData.lastAct-1, :,2, 2,] .= 1
    scenLData.ScenarioLattice.envProbability[envData.firstAct:envData.lastAct-1, :, 1, 2] .=
        tracker[envData.firstAct-1:envData.lastAct-2, :]
    scenLData.ScenarioLattice.envProbability[envData.firstAct:envData.lastAct-1, :, 1, 1] .=
        1 .- tracker[envData.firstAct-1:envData.lastAct-2, :]
    scenLData.ScenarioLattice.envProbability[envData.lastAct, :, :, 1] .= 1
    expandedStates = []
    expandedProb = []
    for t = 1:1:NStage
        if t >= envData.firstAct && t <= envData.lastAct
            newProb = zeros(NStates * 2, NStates * 2)
            if t == envData.lastAct
                newStates = hcat(
                    scenLData.ScenarioLattice.states[t][:, :],
                    repeat([false], NStates),
                )
                newProb[1:NStates, 1:NStates] =
                    scenLData.ScenarioLattice.probability[t][:, :]
                newProb[NStates+1:end, 1:NStates] =
                    scenLData.ScenarioLattice.probability[t][:, :]
            else
                newStates = vcat(
                    hcat(
                        scenLData.ScenarioLattice.states[t][:, :],
                        repeat([false], NStates),
                    ),
                    hcat(
                        scenLData.ScenarioLattice.states[t][:, :],
                        repeat([true], NStates),
                    ),
                )
                for fromState = 1:NStates
                    newProb[fromState, 1:NStates] =
                        scenLData.ScenarioLattice.probability[t][fromState, :] .*
                        (1 - tracker[t-1, fromState])
                    newProb[fromState, NStates+1:end] =
                        scenLData.ScenarioLattice.probability[t][fromState, :] .*
                        tracker[t-1, fromState]
                end
                if t >= envData.firstAct + 1
                    newProb[NStates+1:end, NStates+1:end] =
                        scenLData.ScenarioLattice.probability[t][:, :]
                    #newProb[NStates+1:end,1:NStates] .= 0
                end
            end
        else
            newStates =
                hcat(scenLData.ScenarioLattice.states[t][:, :], repeat([false], NStates))
            newProb = scenLData.ScenarioLattice.probability[t][:, :]
        end
        append!(expandedStates, [newStates])
        append!(expandedProb, [newProb])
    end

    newLattice =
        lattice(expandedStates, expandedProb, scenLData.ScenarioLattice.envProbability)
    return newLattice
end

# -> update probabilities
function Calculate_envConst_prob(envData, data::samplingData, HY::HydroData)

    inflow_lim = envData.actLevel / HY.Scale[envData.envMod]
    t_start = envData.firstAct - 1
    t_end = envData.lastAct - 2

    temp = zeros(size(data.trajectories[:, 1, :]))
    temp[data.trajectories[:, 1, :].>=inflow_lim] .= 1

    count = zeros(size(data.states)[1], size(data.states[1])[1])

    for scen = 1:size(temp)[2]
        for week = t_start:t_end
            if temp[week, scen] > 0
                count[week, Int(data.KmeansClusters[week].assignments[scen])] +=
                    temp[week, scen]
                break
            end
        end
    end

    for week = t_start:t_end
        count[week, :] .= count[week, :] ./ data.KmeansClusters[week].counts
    end

    return count
end

# scale mean price (change average yearly power price profile)
function scalePrice(filename,scenLData)
    meanPrice = Matrix{Float64}(
        CSV.read(joinpath(DataPath, filename), DataFrame, header = false),
    ) #EUR/MWh
    
    #newPricedata =  NormData(scenLData.dataP.data, meanPrice, scenLData.dataP.std);
    newTrajectories = zeros(size(scenLData.trajectories));

    for s = 1:NSamples
        for t= 1:size(scenLData.trajectories)[1]
        newTrajectories[t,1,s] =  scenLData.trajectories[t,1,s]
        newTrajectories[t,2,s] = scenLData.trajectories[t,2,s] .* (meanPrice[t] ./scenLData.dataP.mean[t])
        end
    end

    #states = calculate_states_from_normalized_data(scenLData.KmeansClusters, scenLData.dataI, newPricedata);
    
    temp=[]
    for t = 1:size(scenLData.states)[1]
        updated = zeros(size(scenLData.states[t])) 
        updated[:,1] = scenLData.states[t][:,1]
        updated[:,2] = scenLData.states[t][:,2]  ./ scenLData.dataP.mean[t] .* meanPrice[t]
        append!(temp, [updated])
    end
    
    states = temp
    
    envProbability = zeros(NStage, NStates, 2, 2);
    envProbability[:, :, 1, 1] .= 1;
    newScenarioLattice = lattice(states, scenLData.transitionProb, envProbability);

    newscenLData = samplingData(newTrajectories,
        scenLData.trajectories_normalized,
        scenLData.dataI,
        scenLData.dataP,
        scenLData.KmeansClusters,
        scenLData.transitions,
        scenLData.startProbability,
        scenLData.transitionProb,
        states,
        newScenarioLattice);
    return newscenLData
end


#=
function countProb(T, N, NStates, Inflow, Price, bounds, NStatesPrice)
    probability = zeros(T, NStates, NStates)
    for t = 2:T
        for n = 1:N
            from = [Inflow[t-1, n] Price[t-1, n]]
            to = [Inflow[t, n] Price[t, n]]
            for fromState = 1:NStates
                lower = [bounds[t-1, fromState, 1] bounds[t-1, fromState, 3]]
                upper = [bounds[t-1, fromState, 2] bounds[t-1, fromState, 4]]
                if all(from .>= lower) && all(from .< upper)
                    for toState = 1:NStates
                        lower = [bounds[t, toState, 1] bounds[t, toState, 3]]
                        upper = [bounds[t, toState, 2] bounds[t, toState, 4]]
                        if all(to .>= lower) && all(to .< upper)
                            probability[t, fromState, toState] += 1
                        end
                    end
                end
            end
        end
    end
    return probability
end
=#

#=
function make_lattice_price_inflow(inflow, price, NStates=5, NSamples = 10000)

    #Sample trajectories
    trajectories, trajectories_normalized, dataI, dataP = make_trajectories(inflow, price, NSamples,"Normal")

    NStages = size(trajectories)[1]
    N = size(trajectories)[3]

    #R = [kmeans(trajectories[t,:,:], NStates;maxiter=200) for t =1:T]
    KmeansClusters = [kmeans(trajectories_normalized[t,:,:], NStates;maxiter=200) for t =1:NStages]
    @assert(length(KmeansClusters[1].assignments) == N, "All trajectories are not assigned!")

    transitions = count_transitions(KmeansClusters)
    transitionProb = calculate_probabilities(KmeansClusters, transitions)
    states = calculate_states_from_normalized_data(KmeansClusters, dataI, dataP)

    return lattice(states, transitionProb)
end
=#