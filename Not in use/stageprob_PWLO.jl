#include("interpolate.jl")

struct StageProblem
    model::Any
    res::Any
    spill::Any
    prod::Any
    tank::Any
    disSeg::Any
    inf::Any
    resbalInit::Any
    resbalInitRHS::Any
    resbalStep::Any
    prodeff::Any
    alpha::Any
    gamma::Any
end

function BuildStageProblem(Price)
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
        0.0 <= res[iMod = 1:HY.NMod, iStep = 1:NStep] <= HY.MaxRes[iMod],
        base_name = "res"
    )
    @variable(M, 0.0 <= spi[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "spi")
    @variable(M, 0.0 <= prod[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "prod")
    @variable(M, 0.0 <= tank[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "tank")
    @variable(
        M,
        0.0 <=
        disSeg[iMod = 1:HY.NMod, iSeg = 1:HY.NSeg[iMod], iStep = 1:NStep] <=
        HY.DisMaxSeg[iMod, iSeg],
        base_name = "dseg"
    )
    @variable(M, -AlphaMax <= alpha <= AlphaMax, base_name = "alp")
    @variable(M, resbalInitRHS[iMod = 1:HY.NMod], base_name = "resbalInitRHS")
    @variable(M, inf[iMod = 1:HY.NMod], base_name = "inf")
    #@variable(M,  Alphatab[n=1:NSeg] , base_name="alphaTab")

    @objective(
        M,
        MathOptInterface.MAX_SENSE,
        sum(
            Price[iStep] * prod[iMod, iStep] * NHoursStep - PenTank * tank[iMod, iStep] -
            PenSpi * spi[iMod, iStep] for iMod = 1:HY.NMod for iStep = 1:NStep
        ) + alpha
    )

    #Reservoir balances for first step
    @constraint(
        M,
        resbalInit[iMod = 1:HY.NMod],
        res[iMod, 1] +
        MM3Step * sum(disSeg[iMod, iSeg, 1] for iSeg = 1:HY.NSeg[iMod]) +
        MM3Step * spi[iMod, 1] -
        MM3Step *
        sum(disSeg[uMod, iSeg, 1] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NSeg[uMod]) -
        MM3Step * sum(spi[uMod, 1] for uMod = 1:HY.NUp[iMod]) - MM3Step * tank[iMod, 1] -
        StepFranc * inf[iMod] == resbalInitRHS[iMod]
    )

    #Reservoir balance within stage
    @constraint(
        M,
        resbalStep[iMod = 1:HY.NMod, iStep = 2:NStep],
        res[iMod, iStep] - res[iMod, iStep-1] +
        MM3Step * sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NSeg[iMod]) +
        MM3Step * spi[iMod, iStep] -
        MM3Step * sum(
            disSeg[uMod, iSeg, iStep] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NSeg[uMod]
        ) - MM3Step * sum(spi[uMod, iStep] for umod = 1:Hy.NUp[iMod]) -
        MM3Step * tank[iMod, iStep] - StepFranc * inf[iMod] == 0
    )

    #Hydropower generation
    @constraint(
        M,
        prodeff[iMod = 1:HY.NMod, iStep = 1:NStep],
        prod[iMod, iStep] ==
        sum(HY.Eff[iMod, iSeg] * disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NSeg[iMod])
    )


    #End reservoir in last time step, using SOS2 weights
    #Her er det også mulig å bruke gamma og iterere over hele tabellen

    #Binding the weights to SOS2 and the sum of the weights to be 1

    #use interpolation function from interpolate.jl
    #@expression(M, alpha_bound, linear_iterpolate(res[1,end], Alphatab))
    #@constraint(M, FutureValue, alpha - alpha_bound <= 0 )

    return StageProblem(
        M,
        res,
        spi,
        prod,
        tank,
        disSeg,
        inf,
        resbalInit,
        resbalInitRHS,
        resbalStep,
        prodeff,
        alpha,
        gamma,
    )
end
