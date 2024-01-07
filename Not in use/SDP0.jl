
#SDP algorithm
#Solves using stageproblem.jl
#using Interpolations

struct FutureValueSDP
    ResSeg::Any
    WVTable::Any
    AlphaTable::Any
    DValTable::Any
end

#=
@with_kw struct ResultsSP
    obj
    inflow
    resInit
    alphanfrom
    Reservoir
    Spillage
    Production
    Tanking
    disSeg
    gamma
    totDischarge
end
=#


#AllSeg, ResEnd = twoResStateMatrix(NSeg, HY)
function SDP(DebugSP = false)

    WVTable = zeros(NStage + 1, NScen, NSeg - 1)
    AlphaTable = zeros(NStage + 1, NScen, NSeg)
    DValTable = zeros(NStage + 1, NScen)

    WVTable = zeros(NStage + 1, NScen, NSeg - 1)
    alpha = zeros(MaxIt, NStage, NScen, NSeg)

    #SDP algorithm

    stageprobResult = make_stageprobResult()

    n = 1
    while n <= MaxIt
        println("Iteration: ", n)
        for t = NStage:-1:1
            println("Stage: ", t)

            #if t==NStage
            #    @constraint(SP.model,EndCut[iCut=1:EV.NEndCut],SP.alpha-sum(EV.EndCutCoef[iCut]*SP.res[end]) <= EV.EndCutRHS[iCut])
            #end

            for nfrom = 1:NSeg
                #println("Segment: ", nfrom)
                Alpha = zeros(NScen)
                for scen = 1:NScen
                    #println("Scen: ", i)
                    SP = BuildStageProblem(Price[t, 1:NStep])
                    for iMod = 1:HY.NMod
                        JuMP.fix(SP.inf[iMod], InflowState[iMod, t, scen])
                        JuMP.fix(SP.resbalInitRHS[iMod], ResSeg[nfrom])
                    end

                    #@constraint(SP.model, SP.alpha <= sum(SP.gamma[iSeg]*(WVTable[t+1, i, iSeg-1] * (ResSeg[nfrom]-ResSeg[iSeg-1])+ DValTable[t+1, i]) for iSeg=2:NSeg))
                    @constraint(
                        SP.model,
                        SP.alpha ==
                        sum(SP.gamma[iSeg] * AlphaTable[t+1, scen, iSeg] for iSeg = 1:NSeg)
                    )
                    #@constraint(SP.model, SP.res[1, end] == sum(SP.gamma[iSeg]*ResSeg[iSeg] for iSeg=1:NSeg))

                    #println("Problem nr: ", (NStage-t)*NSeg*NScen + (nfrom-1)*NScen + i )
                    optimize!(SP.model)
                    Alpha[scen] = JuMP.objective_value(SP.model)

                    if DebugSP
                        stageprobResult =
                            check_results(n, t, nfrom, scen, stageprobResult, SP)
                    end


                end
                for fromScen = 1:NScen
                    AlphaTable[t, fromScen, nfrom] = sum(
                        InflowScenarios.InfProb[fromScen, toScen] * Alpha[toScen] for
                        toScen = 1:NScen
                    )
                    if nfrom == 1
                        DValTable[t, fromScen] = AlphaTable[t, fromScen, nfrom]
                    else
                        WVTable[t, fromScen, nfrom-1] =
                            (
                                AlphaTable[t, fromScen, nfrom] -
                                AlphaTable[t, fromScen, nfrom-1]
                            ) / (ResSeg[nfrom] - ResSeg[nfrom-1])
                    end
                end
            end
        end

        diff_matrix = check_diff(WVTable)
        diff = sum(abs.(diff_matrix))
        if diff < conv
            println("Converged iteration:", n, " Diff: ", diff)
            n = MaxIt + 1
        else
            println("Not converged, diff: ", diff)
            WVTable[end, :, :] = WVTable[1, :, :]
            #=
            AlphaTable[end,:,1] = DValTable[1, :]
            AlphaTable[end,:,2:end] = AlphaTable[1, :, 2:end]
            =#

            #Alternativ formulering
            AlphaTable[end, :, 1] .= 0
            for iseg = 1:NSeg-1
                AlphaTable[end, :, iseg+1] =
                    (ResSeg[end] - ResSeg[end-1]) .* WVTable[1, :, iseg] +
                    AlphaTable[end, :, iseg]
            end

            n += 1
        end
    end
    return (FutureValueSDP(ResSeg, WVTable, AlphaTable, DValTable), stageprobResult)
end
