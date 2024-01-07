
#SDP algorithm
#Solves using stageproblem.jl
#using Interpolations

struct FutureValueSDP
    ResSeg::Any
    WVTable::Any
    AlphaTable::Any
    #DValTable
end

function update_state_dependent_input(
    t,
    nfrom,
    iState,
    SP,
    ResSeg,
    ScenarioLattice,
    StepFranc,
    PriceScale,
    NHoursStep,
    HY,
    NStep,
    qMin,
)
    Price = ScenarioLattice.states[t][iState, 2] .* PriceScale
    for iMod = 1:HY.NMod
        JuMP.set_normalized_rhs(
            SP.resbalInit[iMod],
            ResSeg[nfrom][iMod] +
            StepFranc * ScenarioLattice.states[t][iState, 1] * HY.Scale[iMod],
        )

        for iStep = 1:NStep
            set_objective_coefficient(
                SP.model,
                SP.prod[iMod, iStep],
                NHoursStep * Price[iStep],
            )
            
            if !(qMinDependent==0)
                JuMP.set_normalized_rhs(
                    SP.q_min[iMod, iStep], qMin[iMod],
                )
            end
          
            if iStep > 1
                JuMP.set_normalized_rhs(
                    SP.resbalStep[iMod, iStep],
                    StepFranc * ScenarioLattice.states[t][iState, 1] * HY.Scale[iMod],
                )
            end
        end
    end
    return SP
end

function update_future_value_constraint(SP, AlphaTable, NSeg, t, iState)

    if HY.NMod == 2
        for iSegL = 1:NSeg
            for iSegU = 1:NSeg
                JuMP.set_normalized_coefficient(
                    SP.AlphaCon,
                    SP.gamma[iSegU, iSegL],
                    -AlphaTable[t+1, iState, (NSeg*(iSegL-1)+iSegU)],
                )
            end
        end
    elseif HY.NMod == 1
        for iSeg = 1:NSeg
            JuMP.set_normalized_coefficient(
                SP.AlphaCon,
                SP.gamma[iSeg],
                -AlphaTable[t+1, iState, iSeg],
            )
        end
    else
        error("Number of modules not supported.")
    end
    return SP
end

#"SDP algorithm, two dim"
function SDP(
    scenarioLatticeData,
    ResSeg,
    envConst,
    warmstart = 0,
    solveMIP = true,
    useDiag =false,
    DebugSP = false, 
)

    if envConst
        NEnvStates = 2
    else
        NEnvStates = 1
    end

    ScenarioLattice = scenarioLatticeData.ScenarioLattice

    AlphaTable = zeros(NStage + 1, NStates * NEnvStates, length(ResSeg))
    WVTable = zeros(NStage + 1, NStates * NEnvStates, length(ResSeg), HY.NMod)
 
    if warmstart != 0
        println("Using previous watervalues as startpoint..")
        if size(warmstart)[2] == size(WVTable)[2]
            WVTable = warmstart
            WVTable[end,:,:,:,:] = warmstart[1,:,:,:,:]
            AlphaTable[end,:,:,:] = update_futureValue(WVTable, ResSeg, NEnvStates)
         
        else
            WVTable[:,1:NStates,:,:,:] = warmstart
            WVTable[end,1:NStates,:,:,:] = warmstart[1,:,:,:,:]
            AlphaTable[end,:,:,:] = update_futureValue(WVTable, ResSeg, NEnvStates)
           
        end
    end

    #SDP algorithm
    if DebugSP
        stageprobResult = make_stageprobResult()
    end

    #diff_list= []
    println("\nRunning SDP algoritm.")
    println("Number of stages: ", NStage)
    println("Number of scenario states per week: ", NStates)
    println("Number of environmental states per week: ", NEnvStates)
    println("Number of reservoirs: ", HY.NMod)
    println("Number of segments per reservoir: ", NSeg, " total: ", NSeg^HY.NMod)
    println(
        "Number of problems to solve per stage: ",
        NEnvStates * NStates * (NSeg^HY.NMod),
    )
    println(
        "Min number of problems solved per iteration: ",
        NStage * NStates * (NSeg^HY.NMod),
    )
    println(
        "Max number of problems solved per iteration: ",
        NStage * NEnvStates * NStates * (NSeg^HY.NMod),
    )

    n = 1
    SP = @timeit to "build problems" BuildProblem(HY)

    notConvex = false
    MIP_counter = 0
    MIP_counter_t = 0
    nProblems = 0
    optimal = true
    while n <= MaxIt
        optimal = true
        println("Iteration: ", n)
        @time begin
            for t = NStage:-1:1
                @time begin
                    MIP_counter_t = 0
                    NStates = size(ScenarioLattice.states[t])[1]
                    if t > 1
                        NStatesFrom = size(ScenarioLattice.states[t-1])[1]
                    else
                        NStatesFrom = size(ScenarioLattice.states[end])[1]
                    end

                    for nfrom = 1:length(ResSeg) #reservoir state start of stage
                        Alpha = zeros(NStates)
                        for iState = 1:NStates #inflow in stage
                            
                            qMin = copy(HY.qMin)
                            if !(qMinDependent==0)                  
                                if ScenarioLattice.states[t][iState, 1]* HY.Scale[qMinDependent.qMod]*qMinDependent.flowScale[t] < qMinDependent.flowMin[t]
                                    qMin[qMinDependent.qMod] = qMinDependent.qMin #m3/s
                                end
                            end

                            SP =
                                @timeit to "Update input to SP 1" update_state_dependent_input(
                                    t,
                                    nfrom,
                                    iState,
                                    SP,
                                    ResSeg,
                                    ScenarioLattice,
                                    StepFranc,
                                    PriceScale,
                                    NHoursStep,
                                    HY,
                                    NStep,
                                    qMin,
                                )
                          

                            @timeit to "Update  input to SP 2: " begin
                                notConvex = isNotConvex2(WVTable[t+1, iState, :, :])
                                if solveMIP
                                    if (notConvex && useDiag) || (notConvex && n > 1) || (notConvex && warmstart != 0)
                                        if HY.NMod == 2
                                            SOS2_test = @constraint(
                                                SP.model,
                                                [iMod = 1:HY.NMod],
                                                SP.beta[iMod, :] in
                                                MOI.SOS2(collect(1:NSeg))
                                            )
                                            SOS2_diag = @constraint(
                                                SP.model,
                                                SP.χ[:] in MOI.SOS2(collect(1:(2*NSeg-1)))
                                            )
                                        elseif HY.NMod == 1
                                            SOS2_test = @constraint(
                                                SP.model,
                                                [iMod = 1:HY.NMod],
                                                SP.gamma[:] in MOI.SOS2(collect(1:NSeg))
                                            )
                                        end     
                                        MIP_counter += 1
                                        MIP_counter_t += 1      
                                    elseif useDiag
                                        SOS2_diag = @constraint(
                                            SP.model,
                                            SP.χ[:] in MOI.SOS2(collect(1:(2*NSeg-1)))
                                        )
                                        MIP_counter += 1
                                        MIP_counter_t += 1
                                    end
                                end

                                if envConst
                                    for iMod = 1:HY.NMod
                                        SP = relax_dichargeLimit(SP,iMod)
                                        SP = relax_minRes(SP,iMod)
                                        SP = relax_minRes_init(SP,iMod)
                                    end
                                    for envData in envDataList
                                        if t < envData.lastAct && envData.firstAct > 0 &&
                                        Bool(ScenarioLattice.states[t][iState, 3])
                                            SP = add_dichargeLimit(SP, envData, qMin)
                                            SP = relax_minRes(SP,envData.envMod)
                                            SP = relax_minRes_init(SP,envData.envMod)
                                        elseif (t < envData.lastAct && envData.firstAct == 0) ||
                                            (t < envData.lastAct && !Bool(ScenarioLattice.states[t][iState, 3]))
                                            SP = relax_dichargeLimit(SP,envData.envMod)
                                            SP = relax_minRes(SP,envData.envMod)
                                            SP = relax_minRes_init(SP,envData.envMod)
                                        elseif t >= envData.lastAct &&
                                            t <= envData.lastNoDecrease
                                            SP = activate_EnvConstraint(
                                                SP,
                                                t,
                                                nfrom,
                                                iState,
                                                ScenarioLattice,
                                                HY,
                                                ResSeg,
                                                envData,
                                                qMin,
                                            )
                                        end
                                    end
                                end

                                SP = update_future_value_constraint(
                                    SP,
                                    AlphaTable,
                                    NSeg,
                                    t,
                                    iState,
                                )
                            end #timer

                            nProblems += 1
                            @timeit to "optimising problems" optimize!(SP.model)
                            if termination_status(SP.model) == MOI.TIME_LIMIT
                                if MOI.get(optimizer, MOI.ResultCount()) > 0
                                    println(
                                        "Time limit reached. Feasible solution found. t: ",
                                        t,
                                        ", nfrom:",
                                        nfrom,
                                        " , iState",
                                        istate,
                                    )
                                else
                                    println(
                                        "Time limit reached. No feasible solution found. t: ",
                                        t,
                                        ", nfrom:",
                                        nfrom,
                                        " , iState",
                                        istate,
                                    )
                                    optimal = false
                                end
                            elseif termination_status(SP.model) != MOI.OPTIMAL
                                println("NOT OPTIMAL: ", termination_status(SP.model))
                                optimal = false
                            end
                            Alpha[iState] = JuMP.objective_value(SP.model)

                            if DebugSP
                                stageprobResult = check_results(
                                    1,
                                    t,
                                    nfrom,
                                    iState,
                                    stageprobResult,
                                    HY.NMod,
                                    SP,
                                )
                                #TODO add envstate to results!
                            end

                            @timeit to "Clean-up env. state:" begin
                                if solveMIP && ((notConvex && useDiag) || (notConvex && n > 1)  || (notConvex && warmstart != 0)) 
                                    for iMod = 1:HY.NMod
                                        delete(SP.model, SOS2_test[iMod])
                                        if iMod == 2
                                            delete(SP.model, SOS2_diag)
                                        end
                                    end
                                end

                                if useDiag && solveMIP
                                    if HY.NMod == 2
                                        delete(SP.model, SOS2_diag)
                                    end
                                end
                            end #timer
                        end #iState

                        for fromState = 1:NStatesFrom #inflow,price, env state in t-1
                            AlphaTable[t, fromState, nfrom] =
                                @timeit to "calculate expected future value:" sum([
                                    ScenarioLattice.probability[t][fromState, toState] *
                                    Alpha[toState] for toState = 1:NStates
                                ])

                            WVTable[t, fromState, nfrom, :] =
                                @timeit to "calculate water values" calculateWatervalues(
                                    nfrom,
                                    ResSeg,
                                    AlphaTable[t, fromState, :],
                                )
                        end #fromState
                    end #nfrom
                    if solveMIP
                        println("Stage: ", t, " MIP solved: ", MIP_counter_t)
                    else
                        println("Stage: ", t, " MIP solved as LP: ", MIP_counter_t)
                    end
                end #timer
            end #t

            diff = check_diff(WVTable)

            if all(x -> x < conv, diff)
                println("Converged iteration: ", n, ", max diff: ", max(diff...))
                n = MaxIt + 1
            else
                println("Not converged, max diff: ", max(diff...))
                WVTable[end, :, :, :] = WVTable[1, :, :, :]
                AlphaTable[end, :, :] = update_futureValue(WVTable, ResSeg, NEnvStates)
                n += 1

                save(
                    joinpath(Resultspath, string("backup_sdpRes_" * runName, ".jld")),
                    "sdpRes",
                    FutureValueSDP(ResSeg, WVTable, AlphaTable),)
            end
        end #timer
    end #while
    println(MIP_counter, " of ", nProblems, " number of solved problems solved as MIP")
    println("Optimization status: ", optimal)
    if DebugSP
        #return  (FutureValueSDP(ResSeg, WVTable, AlphaTable, DValTable), stageprobResult)
        return (FutureValueSDP(ResSeg, WVTable, AlphaTable), stageprobResult)
    else
        #return  (FutureValueSDP(ResSeg, WVTable, AlphaTable, DValTable))
        return (FutureValueSDP(ResSeg, WVTable, AlphaTable), 0)
    end

end
