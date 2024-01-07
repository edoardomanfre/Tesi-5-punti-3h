function paraSim(
  InputParameters::InputParam,
  ValueTableSDP::FutureValueSDP,
  SimScen::SimScenData,
  runMode,
)
  println("Simulating parallell, ", NSimScen, " scenarios...")

  @unpack ResSeg, WVTable, AlphaTable = ValueTableSDP
  @unpack (scenarios, scenStates) = SimScen
  @unpack NStep, NStage, NSimScen, StepFranc, NHoursStep = InputParameters
  @unpack envConst, solveMIP = runMode

  Eprofit = zeros(NStage)
  Reservoir = zeros(HY.NMod, NSimScen, NStage, NStep)
  Spillage = zeros(HY.NMod, NSimScen, NStage, NStep)
  Production = zeros(HY.NMod, NSimScen, NStage, NStep)
  Q_slack = zeros(HY.NMod, NSimScen, NStage, NStep)
  Res_slack = zeros(HY.NMod, NSimScen, NStage, NStep)
  totDischarge = zeros(HY.NMod, NSimScen, NStage, NStep)
  resInit = zeros(HY.NMod, NSimScen, NStage)
  inflow = zeros(HY.NMod, NSimScen, NStage)
  price = zeros(HY.NMod, NSimScen, NStage, NStep)
  obj = zeros(NSimScen)
  disSeg = []

  if HY.NMod == 1
    gamma = zeros(NSimScen, NStage, NSeg)
  elseif HY.NMod == 2
    gamma = zeros(NSimScen, NStage, NSeg, NSeg)
  end
  for iMod = 1:HY.NMod
    append!(disSeg, [zeros(NSimScen, NStage, NStep, HY.NDSeg[iMod])])
  end

  MIP_counter = 0
  nProblems = 0
  for iScen = 1:NSimScen
    earlyActive_maxDischarge = false
    SP = BuildProblem(InputParameters, HY)
    for t = 1:NStage

      Price = scenarios[iScen][t, 2] .* PriceScale[1:NStep]
      for iMod = 1:HY.NMod
        for iStep = 1:NStep
          set_objective_coefficient(
            SP.model,
            SP.prod[iMod, iStep],
            NHoursStep * Price[iStep],
          )
          if iStep > 1
            JuMP.set_normalized_rhs(
              SP.resbalStep[iMod, iStep],
              StepFranc * scenarios[iScen][t, 1] * HY.Scale[iMod],
            ) #StepFranc[t,1:NStep]
          end
        end


        if t == 1
          JuMP.set_normalized_rhs(
            SP.resbalInit[iMod],
            HY.ResInit0[iMod] + StepFranc * scenarios[iScen][t, 1] * HY.Scale[iMod],
          )
        else
          JuMP.set_normalized_rhs(
            SP.resbalInit[iMod],
            Reservoir[iMod, iScen, t-1, end] +
            StepFranc * scenarios[iScen][t, 1] * HY.Scale[iMod],
          ) 
        end #StepFranc[t,1:NStep]
      end


      if envConst && t >= envDataList[1].firstAct && t < envDataList[1].lastAct
        EnvState = earlyActive_maxDischarge
      else
        EnvState = false
      end

      if envConst
        for envData in envDataList
          SP, earlyActive_maxDischarge = activate_EnvConstraint_sim(
            SP,
            t,
            iScen,
            scenarios,
            HY,
            Reservoir,
            earlyActive_maxDischarge,
            envData,
            qMin,
            NStep,
          )

        end
      end



      SP = update_future_value_constraint_sim(
        SP,
        AlphaTable,
        NSeg,
        t,
        scenStates[iScen][t],
        EnvState,
      )
      #NB: burde jeg interpolere mellom vannverdimatrisene i simuleringen?

      notConvex = isNotConvex(WVTable[t, scenStates[iScen][t], :, :])
      if notConvex
        if solveMIP
          if HY.NMod == 2
            SOS2_test = @constraint(
              SP.model,
              [iMod = 1:HY.NMod],
              SP.beta[iMod, :] in MOI.SOS2(collect(1:NSeg))
            )
            SOS2_diag = @constraint(SP.model, SP.χ[:] in MOI.SOS2(collect(1:(2*NSeg-1))))
          elseif HY.NMod == 1
            SOS2_test = @constraint(
              SP.model,
              [iMod = 1:HY.NMod],
              SP.gamma[:] in MOI.SOS2(collect(1:NSeg))
            )
          end
        end
        MIP_counter += 1
      end

      optimize!(SP.model)
      nProblems += 1
      if termination_status(SP.model) != MOI.OPTIMAL
        println("NOT OPTIMAL: ", termination_status(SP.model))
      end

      obj[iScen] = JuMP.objective_value(SP.model)

      for iMod = 1:HY.NMod
        price[iMod, iScen, t, :] = Price
        for iStep = 1:NStep
          Reservoir[iMod, iScen, t, iStep] = JuMP.value(SP.res[iMod, iStep])
          Spillage[iMod, iScen, t, iStep] = JuMP.value(SP.spill[iMod, iStep])
          Production[iMod, iScen, t, iStep] = JuMP.value(SP.prod[iMod, iStep])
          Q_slack[iMod, iScen, t, iStep] = JuMP.value(SP.q_slack[iMod, iStep])
          Res_slack[iMod, iScen, t, iStep] = JuMP.value(SP.res_slack[iMod, iStep])
          for iSeg = 1:HY.NDSeg[iMod]
            disSeg[iMod][iScen, t, iStep, iSeg] = JuMP.value(SP.disSeg[iMod, iSeg, iStep])
          end
          totDischarge[iMod, iScen, t, iStep] = sum(disSeg[iMod][iScen, t, iStep, :])
        end
      end

      for nU = 1:NSeg
        if HY.NMod == 1
          gamma[iScen, t, nU] = JuMP.value(SP.gamma[nU])
        elseif HY.NMod == 2
          for nL = 1:NSeg
            gamma[iScen, t, nU, nL] = JuMP.value(SP.gamma[nU, nL])
          end
        end
      end
    end
  end
  println("Sim finished")
  println(MIP_counter, " of ", nProblems, " number of solved problems solved as MIP")
  return Results(
    Eprofit,
    Reservoir,
    Spillage,
    Production,
    Q_slack,
    Res_slack,
    disSeg,
    totDischarge,
    resInit,
    inflow,
    price,
    obj,
    gamma,
  )
end