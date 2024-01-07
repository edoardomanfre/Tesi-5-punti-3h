using Distributions
using StatsPlots
using CSV
using Printf
plotly()


struct NormData
    data::Any
    mean::Any
    std::Any
end

#normalize all data
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

function correlation(v1, v2)
    correlation = zeros(size(v1)[2])
    for s = 1:size(v1)[2]
        correlation[s] = cor(v1[:, s], v2[:, s])
    end
    return mean(correlation)
end

function sampleFromDistribution(data, distribution, nSamples)
    if distribution == "Normal"
        dist = MvNormal(vec(mean(data, dims = 1)), cov(data, dims = 1))
    elseif distribution == "LogNormal"
        dist = MvLogNormal(vec(mean(data, dims = 1)), cov(data, dims = 1))
    else
        error("Did not match distribution: use \"Normal\" or \"LogNormal\"")
    end
    samples = rand(dist, nSamples)
    return samples
end

function replaceNegatives(
    data_norm,
    distribution,
    trajectories,
    trajectories_corrected,
    dataI,
    dataP,
    corr_I,
    corr_IP,
)
    negatives = findall(x -> x < 0, trajectories_corrected)

    # Draw new samples for negative values and adjust for correlation, mean and std
    while length(negatives) > 0
        for n in negatives
            s = sampleFromDistribution(data_norm, distribution, 1)
            if n[1] > 1
                trajectories[n[1], 1, n[3]] = s[1] + trajectories[n[1]-1, 1, n[3]] * corr_I
                trajectories[n[1], 2, n[3]] =
                    s[2] + (trajectories[n[1], 1, n[3]] - s[1]) * corr_IP
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

    return trajectories_corrected
end

function makeTrajectories(inflow, price, NSamples = 10000, distribution = "Normal")
    # calculte correlation
    corr_IP = correlation(inflow, price) #average over all scen of the correlation between the inflow and price
    corr_I = correlation(inflow, vcat(inflow[2:end, :], inflow[1, :]')) #average over all scen of the correlation between the inflow in t and t+1

    #Normalize data and calculate parameters
    dataI = normalize(inflow)
    dataP = normalize(price)
    data_norm =
        [reshape(dataI.data, (length(dataI.data))) reshape(dataP.data, (length(dataP.data)))]

    NVariables = size(data_norm)[2]
    trajectories = zeros(size(dataI.data)[1], NVariables, NSamples)
    trajectories_corrected = zeros(size(dataI.data)[1], NVariables, NSamples)

    for t = 1:size(dataI.data)[1]
        trajectories[t, :, :] = sampleFromDistribution(data_norm, distribution, NSamples)
        if t > 1
            temp = trajectories[t, 1, :]
            trajectories[t, 1, :] =
                trajectories[t, 1, :] .+ trajectories[t-1, 1, :] .* corr_I
            trajectories[t, 2, :] =
                trajectories[t, 2, :] .+ (trajectories[t, 1, :] .- temp) .* corr_IP
        end
        trajectories_corrected[t, 1, :] =
            (trajectories[t, 1, :] .* dataI.std[t]) .+ dataI.mean[t]
        trajectories_corrected[t, 2, :] =
            (trajectories[t, 2, :] .* dataP.std[t]) .+ dataP.mean[t]
    end
    trajectories_corrected = replaceNegatives(
        data_norm,
        distribution,
        trajectories,
        trajectories_corrected,
        dataI,
        dataP,
        corr_I,
        corr_IP,
    )


    #print some key parameters desctibing the input data and the sapled data
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
    @printf("correlation inflow: \t\tinput: \t%.2f samples: \t%.2f", corr_I, corr_I_sim)

    return trajectories_corrected
end

struct lattice
    states::Any
    probability::Any
end

function calculateInitialBounds(data, NStates)
    T = size(data)[1]
    Bounds = zeros(T, NStates + 1)
    for t = 1:T
        max = findmax(data[t, :])[1]
        min = findmin(data[t, :])[1]
        Bounds[t, :] = [min + i * (max - min) / (NStates) for i = 0:NStates]
    end
    return Bounds
end

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

function updatePriceBounds(t, bounds, loc, NStates, NStatesPrice)
    locinflow = ceil(loc / NStates)
    locprice = loc - (locinflow - 1) * NStatesPrice
    if locprice == NStatesPrice || (locprice >= NStatesPrice / 2 && (locprice - 2) > 0)
        bounds[t, loc, 3] = bounds[t, loc-1, 3]
        bounds[t, loc-1, 4] = bounds[t, loc-1, 3]
        bounds[t, loc-1, 3] =
            bounds[t, loc-2, 4] - (bounds[t, loc-2, 4] - bounds[t, loc-2, 3]) * 0.25
        bounds[t, loc-2, 4] =
            bounds[t, loc-2, 4] - (bounds[t, loc-2, 4] - bounds[t, loc-2, 3]) * 0.25
    elseif locprice == 1 || (locprice < NStatesPrice / 2 && (locprice + 2) <= NStatesPrice)
        bounds[t, loc, 4] = bounds[t, loc+1, 4]
        bounds[t, loc+1, 3] = bounds[t, loc+1, 4]
        bounds[t, loc+1, 4] =
            bounds[t, loc+2, 3] + (bounds[t, loc+2, 4] - bounds[t, loc+2, 3]) * 0.25
        bounds[t, loc+2, 3] =
            bounds[t, loc+2, 3] + (bounds[t, loc+2, 4] - bounds[t, loc+2, 3]) * 0.25
    else
        bounds[t, loc, 4] =
            bounds[t, loc+1, 3] + (bounds[t, loc+1, 4] - bounds[t, loc+1, 3]) * 0.2
        bounds[t, loc+1, 3] =
            bounds[t, loc+1, 3] + (bounds[t, loc+1, 4] - bounds[t, loc+1, 3]) * 0.2
        bounds[t, loc, 3] =
            bounds[t, loc-1, 4] - (bounds[t, loc-1, 4] - bounds[t, loc-1, 3]) * 0.2
        bounds[t, loc-1, 4] =
            bounds[t, loc-1, 4] - (bounds[t, loc-1, 4] - bounds[t, loc-1, 3]) * 0.2
    end
    return bounds
end

function updateInflowBounds(t, loc, NStates, NStatesPrice, NStatesInflow, bounds)

    locinflow = ceil(loc / NStates)
    locprice = loc - (locinflow - 1) * NStatesPrice
    if locinflow == NStatesInflow || (locinflow >= NStatesInflow / 2 && (locinflow - 2) > 0)
        bounds[t, loc, 3] = bounds[t, loc-NStatesPrice, 3]
        bounds[t, loc-NStatesPrice, 4] = bounds[t, loc-NStatesPrice, 3]
        bounds[t, loc-NStatesPrice, 3] =
            bounds[t, loc-2*NStatesPrice, 4] -
            (bounds[t, loc-2*NStatesPrice, 4] - bounds[t, loc-2*NStatesPrice, 3]) * 0.25
        bounds[t, loc-2*NStatesPrice, 4] =
            bounds[t, loc-2*NStatesPrice, 4] -
            (bounds[t, loc-2*NStatesPrice, 4] - bounds[t, loc-2*NStatesPrice, 3]) * 0.25
    elseif locinflow == 1 ||
           (locinflow < NStatesInflow / 2 && (locinflow + 2) <= NStatesInflow)
        bounds[t, loc, 4] = bounds[t, loc+NStatesPrice, 4]
        bounds[t, loc+NStatesPrice, 3] = bounds[t, loc+NStatesPrice, 4]
        bounds[t, loc+NStatesPrice, 4] =
            bounds[t, loc+2*NStatesPrice, 3] +
            (bounds[t, loc+2*NStatesPrice, 4] - bounds[t, loc+2*NStatesPrice, 3]) * 0.25
        bounds[t, loc+2*NStatesPrice, 3] =
            bounds[t, loc+2*NStatesPrice, 3] +
            (bounds[t, loc+2*NStatesPrice, 4] - bounds[t, loc+2*NStatesPrice, 3]) * 0.25
    else
        bounds[t, loc, 4] =
            bounds[t, loc+NStatesPrice, 3] +
            (bounds[t, loc+NStatesPrice, 4] - bounds[t, loc+NStatesPrice, 3]) * 0.2
        bounds[t, loc+NStatesPrice, 3] =
            bounds[t, loc+NStatesPrice, 3] +
            (bounds[t, loc+NStatesPrice, 4] - bounds[t, loc+NStatesPrice, 3]) * 0.2
        bounds[t, loc, 3] =
            bounds[t, loc-NStatesPrice, 4] -
            (bounds[t, loc-NStatesPrice, 4] - bounds[t, loc-NStatesPrice, 3]) * 0.2
        bounds[t, loc-1, 4] =
            bounds[t, loc-NStatesPrice, 4] -
            (bounds[t, loc-NStatesPrice, 4] - bounds[t, loc-NStatesPrice, 3]) * 0.2
    end
    return bounds
end

function getCombinedBounds(Inflow, Price, NStatesInflow, NStatesPrice)
    bounds = zeros(T, NStates, 4)
    inflowBounds = calculateInitialBounds(Inflow, NStatesInflow)
    priceBounds = calculateInitialBounds(Price, NStatesPrice)
    #add small number to max states
    inflowBounds[:, end] .= inflowBounds[:, end] .+ 0.1
    priceBounds[:, end] .= priceBounds[:, end] .+ 0.1
    #combine into upper and lower boundaries per state variable per state
    for i = 1:NStatesInflow
        for p = 1:NStatesPrice
            bounds[:, (i-1)*NStatesPrice+p, :] .=
                [inflowBounds[:, i] inflowBounds[:, i+1] priceBounds[:, p] priceBounds[
                    :,
                    p+1,
                ]]
        end
    end
    return bounds
end

function makeLattice(NStatesInflow::Int, NStatesPrice::Int, trajectories::Array)

    # Calculate states, combinations of inflow states and price states per week
    for t = 1:T
        inflowStates =
            [(inflowBounds[t, i] + inflowBounds[t, i+1]) / 2 for i = 1:NStatesInflow]
        priceStates = [(priceBounds[t, i] + priceBounds[t, i+1]) / 2 for i = 1:NStatesPrice]

        for i = 1:NStatesInflow
            for p = 1:NStatesPrice
                states[t, (i-1)*NStatesPrice+p, :] = [inflowStates[i] priceStates[p]]
            end
        end
    end

    return lattice(states, probability), bounds
end


#inflow = [1.0 2.1 4.4 1.0 1.4; 4.0 5.2 6.7 3 6; 4.0 5.2 4 6.7 3; 2 1 3 2 2]
#price = [100.0 60.1 10.4 70.0 80.4; 40.0 20.2 10.7 50 20; 30.0 25.2 12 7.7 23; 57 98 65 70 80]

#Read input data
path = "C:\\GitSource\\models\\data\\testmodel"
inflow = Matrix{Float64}(CSV.read(joinpath(path, "inflow.csv"); header = false))
price = Matrix{Float64}(CSV.read(joinpath(path, "price.csv"); header = false))

#Sample trajectories
trajectories = makeTrajectories(inflow, price, 10, "Normal")

T = size(trajectories)[1]
N = size(trajectories)[3]
Inflow = trajectories[:, 1, :]
Price = trajectories[:, 2, :]
NStatesInflow = 3
NStatesPrice = 3
NStates = NStatesInflow * NStatesPrice
probability = zeros(T - 1, NStates, NStates)
states = zeros(T, NStates, 2)

# Calculate bounds for each week
bounds = getCombinedBounds(Inflow, Price, NStatesInflow, NStatesPrice)

# Calculate probabilities
probability = countProb(T, N, NStates, Inflow, Price, bounds, NStatesPrice)

#Adjust bounds
for t = 2:T
    minCount, loc = findmin(sum(probability[t, :, :], dims = 2))
    if minCount < findmin([N / NStates N * 0.1])[1]
        println("Adjusting for week ", t)
        temp = bounds
        bounds = updatePriceBounds(t, bounds, loc[1], NStates, NStatesPrice)
    else
        println("No adjustment week ", t)
    end
    println("Diff week: ", t)
    println(temp[t, :, :] .- bounds[t, :, :])
end




#test, bounds = makeLattice(3, 3, trajectories[1:3,:,:])

#=
# Plot histogram of all drawn samples vs. input data
StatsPlots.histogram(reshape(trajectories[:,1,:], length(trajectories[:,1,:])), normalize=:pdf,  label="samples inflow", color = :blue, alpha=0.3)
StatsPlots.histogram!([inflow...], normalize=:pdf,  label="input inflow", color = :green, alpha=0.3)

StatsPlots.histogram(reshape(trajectories[:,2,:], length(trajectories[:,1,:])), bins=50, normalize=:pdf,  label="samples price", color = :blue, alpha=0.3)
StatsPlots.histogram!([price...], normalize=:pdf,  label="input price", bins=50, color = :green, alpha=0.3)

# Plot sampled data and input for selected week
week = 45
StatsPlots.histogram(reshape(trajectories[week ,1,:], length(trajectories[week,1,:])), bins=20, normalize=:pdf,  label="samples inflow", color = :blue, alpha=0.3)
StatsPlots.histogram!(inflow[week,:], normalize=:pdf, bins=20, label="input inflow", color = :green, alpha=0.3)
StatsPlots.histogram!(reshape(trajectories[week ,2,:], length(trajectories[week,2,:]) ), normalize=:pdf, bins=50, label="samples price", color = :blue, alpha=0.3)
StatsPlots.histogram!(price[week, :], normalize=:pdf, label="input price", color = :green, alpha=0.3)


StatsPlots.plot(trajectories_corrected[:,1,1:100] )
=#

#make multivariat distribution
#Draw samples, correct for mean, std and correlation between periods
#Make lattice with x states
