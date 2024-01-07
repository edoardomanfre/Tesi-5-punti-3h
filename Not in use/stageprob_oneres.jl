#include("interpolate.jl")

struct StageProblem
    model::Any
    res::Any
    spill::Any
    prod::Any
    q_slack::Any
    res_slack::Any
    disSeg::Any
    resbalInit::Any
    resbalStep::Any
    prodeff::Any
    alpha::Any
    gamma::Any
    AlphaCon::Any
    maxRelease::Any
    minReservoirEnd::Any
    minReservoir::Any
    noDecrease_week::Any
end

function BuildStageProblem()
    M = Model(
        with_optimizer(
            CPLEX.Optimizer,
            CPX_PARAM_SCRIND = 0,
            CPX_PARAM_PREIND = 0,
            CPXPARAM_MIP_Tolerances_MIPGap = 0,
        ),
    )

    @variable(
        M,
        0 <= res[iMod = 1:HY.NMod, iStep = 1:NStep] <= HY.MaxRes[iMod],
        base_name = "res"
    )
    @variable(M, 0 <= res_slack[iMod = 1:HY.NMod, iStep = 1:NStep], base_name = "res")
    @variable(M, 0 <= spi[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "spi")
    @variable(M, 0 <= prod[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "prod")
    @variable(M, 0 <= q_slack[iMod = 1:HY.NMod, iStep = 1:NStep], base_name = "q_slack")
    @variable(
        M,
        0 <=
        disSeg[iMod = 1:HY.NMod, iSeg = 1:HY.NDSeg[iMod], iStep = 1:NStep] <=
        HY.DisMaxSeg[iMod, iSeg],
        base_name = "dseg"
    )
    @variable(M, alpha <= Big, base_name = "alp")
    #@variable(M, 0 <= beta[iMod=1:HY.NMod,nU=1:NSeg]  <= 1, base_name="weightReservoir")
    @variable(M, 0 <= gamma[n = 1:NSeg] <= 1, base_name = "weight")

    @objective(
        M,
        MathOptInterface.MAX_SENSE,
        sum(
            prod[iMod, iStep] - 1E10 * q_slack[iMod, iStep] - PenSpi * spi[iMod, iStep] -
            (HY.Eff[iMod, 1] * res_slack[iMod, iStep] / MM3Step) * 39 * 2 for
            iMod = 1:HY.NMod for iStep = 1:NStep
        ) + alpha
    )

    #Reservoir balances for first step
    @constraint(
        M,
        resbalInit[iMod = 1:HY.NMod],
        res[iMod, 1] +
        MM3Step * sum(disSeg[iMod, iSeg, 1] for iSeg = 1:HY.NDSeg[iMod]) +
        MM3Step * spi[iMod, 1] -
        MM3Step *
        sum(disSeg[uMod, iSeg, 1] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NDSeg[uMod]) -
        MM3Step * sum(spi[uMod, 1] for uMod = 1:HY.NUp[iMod]) == 0
    )#resbalInitRHS[iMod] + StepFranc*inf[iMod])

    #Reservoir balance within stage
    @constraint(
        M,
        resbalStep[iMod = 1:HY.NMod, iStep = 2:NStep],
        res[iMod, iStep] - res[iMod, iStep-1] +
        MM3Step * sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NDSeg[iMod]) +
        MM3Step * spi[iMod, iStep] -
        MM3Step * sum(
            disSeg[uMod, iSeg, iStep] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NDSeg[uMod]
        ) - MM3Step * sum(spi[uMod, iStep] for uMod = 1:HY.NUp[iMod]) == 0
    ) #StepFranc*inf[iMod])

    # Lower reservoir volume constraint (slack with punishment)
    @constraint(
        M,
        minResPunish[iMod = 1:HY.NMod, iStep = 1:NStep],
        res[iMod, iStep] + res_slack[iMod, iStep] >= HY.MaxRes[iMod] * 0.1
    )#0.1)

    #Hydropower generation
    @constraint(
        M,
        prodeff[iMod = 1:HY.NMod, iStep = 1:NStep],
        prod[iMod, iStep] ==
        sum(HY.Eff[iMod, iSeg] * disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NDSeg[iMod])
    )

    # Min flow requirment
    @constraint(
        M,
        q_min[iMod = 1:HY.NMod, iStep = 1:NStep],
        sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NDSeg[iMod]) + q_slack[iMod, iStep] >=
        HY.qMin[iMod]
    )


    @constraint(M, AlphaCon, alpha - sum(gamma[iSeg] for iSeg = 1:NSeg) == 0)

    #End reservoir in last time step, using SOS2 weights
    #Her er det også mulig å bruke gamma og iterere over hele tabellen
    @constraint(
        M,
        endResWeight[iMod = 1:HY.NMod],
        res[end] == sum(gamma[iSeg] * ResSeg[iSeg][iMod] for iSeg = 1:NSeg)
    )

    #Binding the weights to SOS2 and the sum of the weights to be 1
    @constraint(M, weightcon, sum(gamma[iSeg] for iSeg = 1:NSeg) == 1)
    #@constraint(M, SOS2[iMod=1:HY.NMod], beta[iMod,:] in MOI.SOS2(collect(1:NSeg)))

    @constraint(
        M,
        maxRelease[iMod = envMod, iStep = 1:NStep],
        sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NDSeg[iMod]) <=
        sum(HY.DisMaxSeg[iMod, iSeg] for iSeg = 1:HY.NDSeg[iMod])
    )

    @constraint(M, minReservoirEnd[iMod = envMod, iStep = NStep], res[iMod, iStep] >= 0)

    @constraint(M, minReservoir[iMod = envMod, iStep = 1:NStep], res[iMod, iStep] >= 0)

    @constraint(M, noDecrease_week[iMod = envMod, iStep = NStep], res[iMod, iStep] >= 0)

    return StageProblem(
        M,
        res,
        spi,
        prod,
        q_slack,
        res_slack,
        disSeg,
        resbalInit,
        resbalStep,
        prodeff,
        alpha,
        gamma,
        AlphaCon,
        maxRelease,
        minReservoirEnd,
        minReservoir,
        noDecrease_week,
    )
end
