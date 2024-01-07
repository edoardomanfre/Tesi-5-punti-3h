
struct DetProblem
  model::Any
  res::Any
  spi::Any
  prod::Any
  disSeg::Any
  β::Any
  resbalInit::Any
  resbalStage::Any
  prodeff::Any
  maxRelease::Any
  minReservoir::Any
  noFlippFloppStage::Any
  noFlippFloppStageFist::Any
  α::Any
end



function BuildDeterministicProblem(HY::HydroData, Price, Inflow, resbalInitRHS)

  #=
  1. No early activation
  2. No "no deacrease"
  3. Only no discharge +  early deactivation

  =#
  M = Model(
    with_optimizer(
      CPLEX.Optimizer,
      CPX_PARAM_SCRIND = 0,
      CPX_PARAM_PREIND = 0,
      CPXPARAM_MIP_Tolerances_MIPGap = 0,
      CPXPARAM_MIP_Tolerances_AbsMIPGap = 0,
    ),
  )

  @variable(M, 0 <= res[iMod = 1:HY.NMod, iStage = 1:NStage] <= HY.MaxRes[iMod], base_name = "res")

  @variable(M, 0 <= spi[iMod = 1:HY.NMod, iStage = 1:NStage] <= Big, base_name = "spi")
  @variable(M, 0 <= prod[iMod = 1:HY.NMod, iStage = 1:NStage] <= Big, base_name = "prod")
  @variable(M, 0 <= q_slack[iMod = 1:HY.NMod, iStage = 1:NStage], base_name = "q_slack")
  @variable(
    M,
    0 <=
    disSeg[iMod = 1:HY.NMod, iSeg = 1:HY.NDSeg[iMod], iStage = 1:NStage] <=
    HY.DisMaxSeg[iMod, iSeg],
    base_name = "dseg"
  )
  @variable(M, β[iMod = 1:HY.NMod, iStage = 1:NStage], Bin, base_name = "alp")
  @variable(M, α[iMod = 1:HY.NMod, iStage = NStage] <= Big, base_name = "alp")
  #@variable(M, 0 <= δ[n=1:NCut]<= 1, base_name="weight")

  @objective(
    M,
    MathOptInterface.MAX_SENSE,
    sum(
      sum(
        NHoursStep * Price[iStage] * prod[iMod, iStage] - PenSpi * spi[iMod, iStage] for
        iStage = 1:NStage
      ) + α[iMod, NStage] for iMod = 1:HY.NMod
    )
  )

  #@constraint(M, AlphaCon[imod=1:HY.NMod, n=1:NCut], α <= cutRHS[n] + δ[n]*res[iMod,NStage])
  @constraint(
    M,
    AlphaCon1[iMod = 1:HY.NMod, iStage = NStage],
    α[iMod, iStage] <= 0 + 55555 * res[iMod, iStage]
  )
  #@constraint(M, AlphaCon0[iMod=1:HY.NMod,iStage=NStage], α[iMod, iStage] <= 7185*10^6*100)
  @constraint(
    M,
    AlphaCon2[iMod = 1:HY.NMod, iStage = NStage],
    α[iMod, iStage] <= (0.3 * HY.MaxRes[iMod]) * (55555 - 7185) + 7185 * res[iMod, iStage]
  )
  @constraint(
    M,
    AlphaCon3[iMod = 1:HY.NMod, iStage = NStage],
    α[iMod, iStage] <=
    (0.3 * HY.MaxRes[iMod]) * (55555 - 7185) +
    0.7 * HY.MaxRes[iMod] * (7185 - 0) +
    0 * res[iMod, iStage]
  )

  #Reservoir balances for first step
  @constraint(
    M,
    resbalInit[iMod = 1:HY.NMod, iStage = 1],
    res[iMod, iStage] +
    MM3Step * sum(disSeg[iMod, iSeg, iStage] for iSeg = 1:HY.NDSeg[iMod]) +
    MM3Step * spi[iMod, iStage] -
    MM3Step *
    sum(disSeg[uMod, iSeg, iStage] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NDSeg[uMod]) -
    MM3Step * sum(spi[uMod, iStage] for uMod = 1:HY.NUp[iMod]) ==
    resbalInitRHS[iMod] + Inflow[iStage] * HY.Scale[iMod]
  )

  #Reservoir balance between stage
  @constraint(
    M,
    resbalStage[iMod = 1:HY.NMod, iStage = 2:NStage],
    res[iMod, iStage] - res[iMod, iStage-1] +
    MM3Step * sum(disSeg[iMod, iSeg, iStage] for iSeg = 1:HY.NDSeg[iMod]) +
    MM3Step * spi[iMod, iStage] -
    MM3Step *
    sum(disSeg[uMod, iSeg, iStage] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NDSeg[uMod]) -
    MM3Step * sum(spi[uMod, iStage] for uMod = 1:HY.NUp[iMod]) ==
    Inflow[iStage] * HY.Scale[iMod]
  )

  #Reservoir balance last stage
  # @constraint(M, resbalStepEnd[iMod=1:HY.NMod, iStage=NStage], res[iMod,iStage]  ==  resbalInitRHS[iMod])

  #Hydropower generation
  @constraint(
    M,
    prodeff[iMod = 1:HY.NMod, iStage = 1:NStage],
    prod[iMod, iStage] ==
    sum(HY.Eff[iMod, iSeg] * disSeg[iMod, iSeg, iStage] for iSeg = 1:HY.NDSeg[iMod])
  )

  # No discharge
  @constraint(
    M,
    maxRelease[iMod = envData.envMod, iStage = envData.lastAct:envData.lastMaxDisch],
    sum(disSeg[iMod, iSeg, iStage] for iSeg = 1:HY.NDSeg[iMod]) <=
    HY.qMin[iMod] +
    (1 - β[iMod, iStage]) * sum(HY.DisMaxSeg[iMod, iSeg] for iSeg = 1:HY.NDSeg[iMod])
  )

  @constraint(
    M,
    minReservoir[iMod = 1:HY.NMod, iStage = envData.lastAct:envData.lastMaxDisch],
    res[iMod, iStage] + β[iMod, iStage] * envData.deactLevel >= envData.deactLevel
  )

  #@constraint(M, minReservoir2[iMod=1:HY.NMod, iStage=envData.lastAct:envData.lastMaxDisch],
  #    res[iMod, envData.lastAct] + β[iMod,iStage]*envData.deactLevel>= envData.deactLevel)

  @constraint(
    M,
    noFlippFloppStage[iMod = 1:HY.NMod, iStage = envData.lastAct+1:envData.lastMaxDisch],
    β[iMod, iStage-1] >= β[iMod, iStage]
  )

  @constraint(
    M,
    noFlippFloppStageFist[iMod = 1:HY.NMod, iStage = envData.lastAct],
    β[iMod, iStage] == 1
  )

  @constraint(
    M,
    beforeConst[iMod = 1:HY.NMod, iStage = 1:envData.lastAct-1],
    β[iMod, iStage] == 0
  )

  @constraint(
    M,
    deactivateConst[iMod = 1:HY.NMod, iStage = envData.lastMaxDisch+1:NStage],
    β[iMod, iStage] == 0
  )

  return DetProblem(
    M,
    res,
    spi,
    prod,
    disSeg,
    β,
    resbalInit,
    resbalStage,
    prodeff,
    maxRelease,
    minReservoir,
    noFlippFloppStage,
    noFlippFloppStageFist,
    α,
  )
end

struct Results
  startResLevels::Any
  Reservoir::Any
  Spillage::Any
  Production::Any
  activeMaxDisch::Any
  disSeg::Any
  totDischarge::Any
  resInit::Any
  obj::Any
  alpha::Any
  #profit::Any
  #weekly_profit::Any
end

function detOpt(SimScen::SimScenData, envConst)
  println("Solving deterministic model...")

  @unpack (scenarios, scenStates) = SimScen
  startResLevels = [i for i = 0:0.1:1]
  Nreslevels = length(startResLevels)

  #profit = zeros(HY.NMod,NSimScen,Nstage)
  #weekly_profit = zeros(Hy.NMod,NSimScen,Nstage)
  Reservoir = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  Spillage = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  Production = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  totDischarge = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  activeMaxDisch = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  resInit = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  inflow = zeros(HY.NMod, NSimScen, Nreslevels, NStage)
  obj = zeros(NSimScen, Nreslevels)
  disSeg = []
  alpha = zeros(HY.NMod, Nreslevels, 1)

  for iMod = 1:HY.NMod
    append!(disSeg, [zeros(NSimScen, Nreslevels, NStage, HY.NDSeg[iMod])])
  end

  for iScen = 1:NSimScen
    not_active_maxDischarge = true
    Price = scenarios[iScen][:, 2]
    Inflow = scenarios[iScen][:, 1]
    for startRes = 1:Nreslevels
      resbalInitRHS = zeros(HY.NMod)
      for iMod = 1:HY.NMod
        resbalInitRHS[iMod] = HY.MaxRes[iMod] * startResLevels[startRes]
      end
      SP = BuildDeterministicProblem(HY, Price, Inflow, resbalInitRHS)
      optimize!(SP.model)
      if termination_status(SP.model) != MOI.OPTIMAL
        println(
          "Scen ",
          iScen,
          ", startRes ",
          startRes,
          ": NOT OPTIMAL, ",
          termination_status(SP.model),
        )
      end

      obj[iScen, startRes] = JuMP.objective_value(SP.model)

      for iMod = 1:HY.NMod
        for iStage = 1:NStage
          Reservoir[iMod, iScen, startRes, iStage] = JuMP.value(SP.res[iMod, iStage])
          Spillage[iMod, iScen, startRes, iStage] = JuMP.value(SP.spi[iMod, iStage])
          Production[iMod, iScen, startRes, iStage] = JuMP.value(SP.prod[iMod, iStage])
          activeMaxDisch[iMod, iScen, startRes, iStage] = JuMP.value(SP.β[iMod, iStage])
          for iSeg = 1:HY.NDSeg[iMod]
            disSeg[iMod][iScen, startRes, iStage, iSeg] =
              JuMP.value(SP.disSeg[iMod, iSeg, iStage])
          end
          totDischarge[iMod, iScen, startRes, iStage] =sum(disSeg[iMod][iScen, startRes, iStage, :])
          #profit[iMod,iScen,t,iStep] = price[iMod,iScen,t,iStep]*Production[iMod,iScen,t,iStep]
          #weekly_profit[iMod,iScen,t] = weekly_profit[iMod,iScen,t]+profit[iMod,iSce,t,iStep]
        end #iStage
        alpha[iMod, startRes, 1] = JuMP.value(SP.α[iMod, NStage])
      end #iMod
    end # startRes
    println("Scen ", iScen, ": optimised.")
  end # SimScen

  println("Optimisation of all instances finished")
  return Results(
    startResLevels,
    Reservoir,
    Spillage,
    Production,
    activeMaxDisch,
    disSeg,
    totDischarge,
    resInit,
    obj,
    alpha,
    #profit,
    #weekly_profit
  )
end
