function make_stageprobResult()
    MaxIt = 1
    stageprobResult = Dict(
        "obj" => zeros(MaxIt, NStage, length(ResSeg), NStates),
        #"inflow" => zeros(MaxIt, NStage, length(ResSeg), NStates, HY.NMod),
        #"resInit" => zeros(MaxIt, NStage, length(ResSeg), NStates, HY.NMod),
        "alphanfrom" => zeros(MaxIt, NStage, length(ResSeg), NStates),
        "Reservoir" => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "Spillage" => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "Production" => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "Q_slack" => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "disSeg" => [
            zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NDSeg[iMod]) for
            iMod = 1:HY.NMod
        ],
        "totDischarge" => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
    )

    if HY.NMod == 1
        stageprobResult["gamma"] = zeros(MaxIt, NStage, length(ResSeg), NStates, NSeg)
    else
        stageprobResult["gamma"] = zeros(MaxIt, NStage, length(ResSeg), NStates, NSeg, NSeg)
    end
    return stageprobResult
end

function check_results(n, t, nfrom, state, stageprobResult, NMod, SP::StageProblem)
    stageprobResult["obj"][n, t, nfrom, state] = JuMP.objective_value(SP.model)
    stageprobResult["alphanfrom"][n, t, nfrom, state] = JuMP.value(SP.alpha)
    for iMod = 1:NMod
        #stageprobResult["inflow"][n, t, nfrom, state, iMod] =
        #    JuMP.value(SP.inf[iMod])
        #stageprobResult["resInit"][n, t, nfrom, state, iMod] =
        #    JuMP.value(SP.resbalInitRHS[iMod])
        for iStep = 1:NStep
            stageprobResult["Reservoir"][n, t, nfrom, state, iStep, iMod] =
                JuMP.value(SP.res[iMod, iStep])
            stageprobResult["Spillage"][n, t, nfrom, state, iStep, iMod] =
                JuMP.value(SP.spill[iMod, iStep])
            stageprobResult["Production"][n, t, nfrom, state, iStep, iMod] =
                JuMP.value(SP.prod[iMod, iStep])
            stageprobResult["Q_slack"][n, t, nfrom, state, iStep, iMod] =
                JuMP.value(SP.q_slack[iMod, iStep])
            for iSeg = 1:HY.NDSeg[1]
                stageprobResult["disSeg"][iMod][n, t, nfrom, state, iStep, iSeg] =
                    JuMP.value(SP.disSeg[iMod, iSeg, iStep])
            end
            stageprobResult["totDischarge"][n, t, nfrom, state, iStep, iMod] =
                sum(stageprobResult["disSeg"][iMod][n, t, nfrom, state, iStep, :])
            # stageprobResults["profit"][n, t, nfrom,state,iStep,iMod]=stageprobResult["Production"][n, t, nfrom, state, iStep, iMod]*stageprobResult["price"][n, t, nfrom, state, iStep, iMod]
            # stageprobResults["weekly_profit"][n, t, nfrom,iStep,iMod]
        end
    end

    for nU = 1:NSeg
        if HY.NMod == 1
            stageprobResult["gamma"][n, t, nfrom, state, nU] = JuMP.value(SP.gamma[nU])
        elseif HY.NMod == 2
            for nL = 1:NSeg
                stageprobResult["gamma"][n, t, nfrom, state, nU, nL] =
                    JuMP.value(SP.gamma[nU, nL])
            end
        end
    end
    return stageprobResult
end
