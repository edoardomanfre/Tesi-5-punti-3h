#---------------------------------------------------------------------------------------------
# Description: Solving the scheduling problem by SDP
#              Repeatedly modify and solves an optimisation problem defined in stageprob.mod
#              Computes water-value tables for each stage and state
#
# Developers: LInn Emelie Schäffer
#
# Assumptions: 
#  -State variables: {reservoir} + {stochastic processes} 
#                    + {option:Envionmental state}
#  -Inflow and price sampled together (to preserve correlation) 
#   and is asjusted for autocorrelation 
# 
# Handling of state-dependent constrains:
#
#
#---------------------------------------------------------------------------------------------

# SDP algorithm, two reservoir system
function SDP(
  InputParameters,
  SolverParameters,       #CPLEX
  HY,                     #Hydraulic data dei due power system
  scenarioLatticeData,    #Tutti i possibili scenari generati con valori di inflow e prezzi
  ResSeg,                 #5x5=25 possibili combinazioni di volumi
  PriceScale,             #Matrice 52x7 o 52x56: 7gg o 8x7 ore
  envDataList,            #Lista dei vincoli ambientali con relativi valori
  runMode,
  FinalResPath,
  warmstart = 0,          #Calcola tutto dall'inizio. Se non 0 prende valori di WV già salvati
 )

  # input parameters
  @unpack (NSeg, NStage, NStates, MaxIt, conv, StepFranc, NHoursStep, NStep, LimitPump) =    InputParameters #StepFlow

  @unpack (envConst, solveMIP, DebugSP, ramping_constraints) = runMode

  if envConst
    NEnvStates = 2
  else
    NEnvStates = 1
  end

  ScenarioLattice = scenarioLatticeData.ScenarioLattice

  # make arrays for results
  AlphaTable = zeros(NStage + 1, NStates * NEnvStates, length(ResSeg))                     # Creo una matrice nulla 53x5(o10)x25
  WVTable = zeros(NStage + 1, NStates * NEnvStates, length(ResSeg), HY.NMod)               # Creo una matrice 53x5(o10)x25x2 (una per upper e l'altra per lower reservoir)

  if warmstart != 0
    println("Using previous watervalues as startpoint..")                                  # Se ho gia' dei valori di WVT posso cominiciare gia' da quelli invece di fare un'altra matrice nulla
    if size(warmstart)[2] == size(WVTable)[2]                                              
      WVTable = warmstart
      WVTable[end, :, :, :, :] = warmstart[1, :, :, :, :]
      AlphaTable[end, :, :, :] = update_futureValue(WVTable, ResSeg, NEnvStates)

    else
      WVTable[:, 1:NStates, :, :, :] = warmstart
      WVTable[end, 1:NStates, :, :, :] = warmstart[1, :, :, :, :]
      AlphaTable[end, :, :, :] = update_futureValue(WVTable, ResSeg, NEnvStates)           # Se il codice crascia, posso iniizare da valori gia' esistenti
    end
  end


  if DebugSP
    # save results before convergence
    stageprobResult = make_stageprobResult() 
  end

  # print information about problem
  println("\nRunning SDP algoritm.")
  println("Number of stages: ", NStage)
  println("Number of scenario states per week: ", NStates)
  println("Number of environmental states per week: ", NEnvStates)
  println("Number of reservoirs: ", HY.NMod)
  #println("Number of segments per reservoir: ", NSeg, " total: ", NSeg^HY.NMod)
  println("Number of segment per upper reservoir:", NSeg[1], " and lower ", NSeg[2], " total:", NSeg[1]*NSeg[2])
  #println("Number of problems to solve per stage: ", NEnvStates * NStates * (NSeg^HY.NMod))
  println("Number of problems to solve per stage: ", NEnvStates * NStates * (NSeg[1]*NSeg[2]))
  #println("Min number of problems solved per iteration: ",NStage * NStates * (NSeg^HY.NMod),)
  println("Min number of problems solved per iteration: ",NStage * NStates * (NSeg[1]*NSeg[2]),)
  #println("Max number of problems solved per iteration: ",NStage * NEnvStates * NStates * (NSeg^HY.NMod),)
  println("Max number of problems solved per iteration: ",NStage * NEnvStates * NStates * (NSeg[1]*NSeg[2]),)

  # START SDP ALGORITHM
  # --------------

  # build initial optimisation problem
  SP = @timeit to "build problems" BuildProblem(InputParameters, HY, SolverParameters)                                        # Setup the model (stageprob) considering 1 or 2 reservoirs, writing the variable,objjective f, constraints etc

  # set counters
  MIP_counter = 0                                                                                                             # Conta il numero delle volte che vado in un modello che devo risolvere
  nProblems = 0
  #MIP_counter_t = 0
  #optimal = true
  #notConvex = false
  n = 1
  optimal = true
 
  while n <= MaxIt                                                          # Comincio ad iterare (MaxIt=100) fino a convergenza o numero massimo di iterazioni                                                                                        
    optimal = true
    println("Iteration: ", n)

    @time begin
      for t = NStage:-1:1                                                   # Itero dalla settimana 52 e vado indietro                                                                                        
        @time begin
          MIP_counter_t = 0                                                 # Quanti MIP sono stati risolti per a settimana

          add_dischargeLimitPump = false                                                                                                 

          NStatesTo = size(ScenarioLattice.states[t])[1]                    # total number of stochastic and environmental states in current stage (5 or 10): 5 possibili stati per ogni settimana                                                                            # I have 5 states for each week - in alcuni casi ho l'espansione 
          
          for nfrom = 1:length(ResSeg)                                      # Itero per tutte le combinazioni di volume
            # for all reservoir states (at beginning of stage)                                                                 
            Alpha = zeros(NStatesTo)                                        # Create the matrix Alpha with dim.5 
            
            for iState = 1:NStatesTo                                        # Comincio un ciclo per tutti gli stati del sistema = 1:5                                                              
              # for all stochastic states (inflow, price, env)        
              @timeit to "Update input to decision problem : " begin

                # PRE-PROCESSING OF DECISION PROBLEM
                # ----------------------------------
                                                                          
                # update optimisation problem                               # Per tutte le settimane, combinazioni di reservoir e stati del sistema aggiorno prezzi e inflow                                                           
                SP = @timeit to "Update decision problem" update_state_dependent_input(       
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
                  HY.qMin,
                  ramping_constraints,
                  #StepFlow,
                )


                # check if future value table is convex or not     #Controlla se la funzione WV (derivata di alphavalue) è convessa. Se non lo è viene linearizzata con SOS2
                @timeit to "Check convex: " notConvex =
                  isNotConvex(WVTable[t+1, iState, :, :])                                                                                 # Function present in "waterValueFunctions"

                @timeit to "Set SOS2 requirment: " begin                                                                                  # allowing problem to be a MIP problem, if false it doesn't continue

                  # adding SOS2 constraint if non-linear and MIP allowed
                  if solveMIP
                    if (notConvex && n > 1) || (notConvex && warmstart != 0)                                                              # if not convex, can solve but not in first iteration - start from second
                      if HY.NMod == 2
                        SOS2_test_upper = @constraint(
                          SP.model,
                         # [iMod = 1],
                          SP.beta_upper[:] in MOI.SOS2(collect(1:NSeg[1]))       # separate for Upper and Lower
                        )
                        SOS2_test_lower = @constraint(
                          SP.model,
                        #  [iMod = 2],
                          SP.beta_lower[:] in MOI.SOS2(collect(1:NSeg[2]))       # separate for Upper and Lower
                        )

                        #  SOS2_diag =                                                                                 #togliere SOS2_diag
                       #   @constraint(SP.model, SP.χ[:] in MOI.SOS2(collect(1:(2*NSeg-1))))
                      
                      elseif HY.NMod == 1
                        SOS2_test = @constraint(
                          SP.model,
                        #  [iMod = 1],
                          SP.gamma[:] in MOI.SOS2(collect(1:NSeg[1])) # NSeg[iMod]
                        )
                      end
                      MIP_counter += 1
                      MIP_counter_t += 1
                    end
                  end
                end #timer SOS2

                @timeit to "Set env requirements: " begin                              #Impongo i vincoli ambientali: state dependent constraint
                  # check environmental requirements and update problem
                  
                  if envConst                                                          #Se ci sono i vincoli a seconda della settimana vengono imposti i limiti di discharge e volume

                    for envData in envDataList
                      if t < envData.lastAct && envData.firstAct > 0 && Bool(ScenarioLattice.states[t][iState, 3])      #Se 0<t<23 
                        SP = add_dichargeLimit(SP, envData, NStep)                              
                        add_dischargeLimitPump = true                                                                   #Impongo limite sul pompaggio (52)
                        #SP = add_dischargeLimitPump(SP,NStep)
                        SP = relax_minRes(SP, envData.envMod, NStep)
                        SP = relax_minRes_init(SP, envData.envMod, NStep)
                      elseif (t < envData.lastAct && envData.firstAct == 0) || (t < envData.lastAct && !Bool(ScenarioLattice.states[t][iState, 3]))
                        SP = relax_dichargeLimit(SP, envData.envMod, NStep)
                        add_dischargeLimitPump = false
                        #SP = relax_dischargeLimitPump(SP,NStep)
                        SP = relax_minRes(SP, envData.envMod, NStep)
                        SP = relax_minRes_init(SP, envData.envMod, NStep)
                      elseif t >= envData.lastAct && t <= envData.lastNoDecrease       #If 23<t<38 , activates the Environmental constraints
                        SP, add_dischargeLimitPump = activate_EnvConstraint(
                          SP,
                          t,
                          nfrom,
                          iState,
                          ScenarioLattice,
                          HY,
                          ResSeg,
                          envData,
                          NStep,
                          add_dischargeLimitPump,
                        )
                      end
                    end

                  end
                end #timer env

                if HY.Station_with_pump==1 && !add_dischargeLimitPump                         # if is false, there are no restric.pumping due to state-constraints on max discharge
                  SP = DeactivationPump_SDP(SP,HY,ResSeg,LimitPump,nfrom,NStep)    
                elseif HY.Station_with_pump==1 && add_dischargeLimitPump                      # true - add the constraint on pumping!
                  SP = add_disLimitPump(SP,NStep)
                end

                @timeit to "Update future value: " SP =
                  update_future_value_constraint(SP, AlphaTable, NSeg, t, iState)
                end #timer

              # SOLVE DECISION problem
              # ----------------------

              nProblems += 1
              @timeit to "optimising problems" optimize!(SP.model)                            # Start optimization of the problem
              if termination_status(SP.model) == MOI.TIME_LIMIT                               # Se il modello e' giunto a termine  
                if MOI.get(optimizer, MOI.ResultCount()) > 0
                  println(
                    "Time limit reached. Feasible solution found. t: ",                       # Se ha trovato una soluzione ottima
                    t,
                    ", nfrom:",
                    nfrom,
                    " , iState",
                    istate,
                  )
                else
                  println(
                    "Time limit reached. No feasible solution found. t: ",                    # Se non ha trovato una soluzione fattibile
                    t,
                    ", nfrom:",
                    nfrom,
                    " , iState",
                    istate,
                  )
                  optimal = false
                end
              elseif termination_status(SP.model) != MOI.OPTIMAL
                println("NOT OPTIMAL: ", termination_status(SP.model))                        # Se non ha trovato alcuna soluzione (alcun valore ottimale)
                optimal = false
              end

              # SAVE RESULTS 
              # ------------

              Alpha[iState] = JuMP.objective_value(SP.model)             #Salvo i risultati per tutti gli stati (5): ottengo i valori di expected future profit
              

              if DebugSP
                stageprobResult =
                  check_results(1, t, nfrom, iState, stageprobResult, HY.NMod, SP)     # Function in "getResults"- save the results for Eprofit, obj , reservoir, spillage, production
                #TODO add envstate to results!
              end

              # CLEAN-UP PROBLEM                                                # Se non rimossi, puo' dare problemi nel aggiungerli - solo per renderlo piu' leggero ed evitare vincoli inutili
              # ----------------                                                # e per preparalo per la prossima iterazione

              @timeit to "Clean-up env. state:" begin                           #Rimuovo tutti i constraints attivati che poi verranno rimessi o no nella prossima iterazione
                # remove SOS2 constraints
                if solveMIP && ((notConvex && n > 1) || (notConvex && warmstart != 0)) #((notConvex && useDiag) || 
                  if HY.NMod==2
                    delete(SP.model, SOS2_test_upper)
                    delete(SP.model, SOS2_test_lower)                    
                  elseif HY.NMod==1
                    delete(SP.model,SOS2_test)
                  end
                end

                # remove environmental constraints and pump constraints                                                                             # Assicurarsi di non portarsi dietro vincoli ambientali inutili
                if envConst
                  for iMod = 1:HY.NMod
                    SP = relax_dichargeLimit(SP, iMod, NStep)
                    #add_dischargeLimitPump = false
                    SP = relax_minRes(SP, iMod, NStep)
                    SP = relax_minRes_init(SP, iMod, NStep)
                    if HY.Station_with_pump==1
                      SP = relax_disLimitPump(SP,NStep)
                    end
                  end
                end

                # remove pumping constraints if there are no env.constraints
                if !envConst
                  if HY.Station_with_pump==1
                    SP = relax_disLimitPump(SP,NStep)
                  end
                end

              end #timer
            end #iState

            # CALCULATE EXPECTED VALUE  
            # ------------------------

            if t > 1
              NStatesFrom = size(ScenarioLattice.states[t-1])[1]
            else
              NStatesFrom = size(ScenarioLattice.states[end])[1]
            end
            
            for fromState = 1:NStatesFrom #inflow,price, env state in t-1
              
              AlphaTable[t, fromState, nfrom] =                            #[52, 5, 25]. Dopo aver calcolato Alpha calcolo AlphaTable=sommatoria probability(t)[from state, to state]+Alpha[to state]                                                 # for week t, State of the system (1-5), for nfrom = combination of reservoir segments
                @timeit to "calculate expected future value:" sum([
                  ScenarioLattice.probability[t][fromState, toState] * Alpha[toState] for
                  toState = 1:NStates
                ])

              WVTable[t, fromState, nfrom, :] =                           #[week, state,res combination, upper or lower]
                @timeit to "calculate water values" calculateWatervalues(                                                   # function in line 458
                  nfrom,
                  ResSeg,
                  AlphaTable[t, fromState, :],                           # Considering Alpha in week t, at state "fromState" and for all resSeg combinations
                )
            end #fromState                                               # Finishes when we have done it for all states in week t
          
          end #nfrom                                                     # Finishes when done for all seg.combinations                        

          if solveMIP
            println("Stage: ", t, " MIP solved: ", MIP_counter_t)        # Quanti MIP ho risolto
            #println("Alpha week 53: ", AlphaTable[end,1,1:5])
            #println("Alpha state $t", AlphaTable[t,1,1:5])
            #println("Difference:", AlphaTable[t+1,1,1:5]-AlphaTable[t,1,1:5])
          else
            println("Stage: ", t, " MIP solved as LP: ", MIP_counter_t)
          end
        end #timer
      end #t

      # CHECK FOR CONVERGENCE 
      # ----------------------

      # Calculate difference in water values
      @timeit to "calculate diff: " diff = check_diff(WVTable)                         # Function in line 444 - returns the difference

      # check towards convergence criteria
      @timeit to "check conv. and update: " begin
        if all(x -> x < conv, diff)                                                    # if diff <= conv (0.01) we have reached convergence
          println("Converged iteration: ", n, ", max diff: ", max(diff...))
          n = MaxIt + 1
        else
          println("Not converged, max diff: ", max(diff...))                           # in the case we didn't reach convergence, we update the WVTable and AlphaTable for next iteration

          # update values for next iteration
          WVTable[end, :, :, :] = WVTable[1, :, :, :]
          AlphaTable[end, :, :] = update_futureValue(WVTable, ResSeg, NEnvStates)
          n += 1

          # store temporary results (in case of code-crash)
          @timeit to "save results iteration: " save(
            joinpath(FinalResPath, string("backup_sdpRes.jld")),
            "sdpRes",
            FutureValueSDP(ResSeg, WVTable, AlphaTable),
          )
        end

      end #timer check convergence
    end #timer
  end #while                                                                                                               # end - start a new iteration if needed

  # STORE FINAL RESULTS FROM SDP ALGORITHM 
  # ---------------------------------------

  # print final optimisation status
  println(MIP_counter, " of ", nProblems, " number of solved problems solved as MIP")
  println("Optimization status: ", optimal)

  # store results to file
  if DebugSP
    return (FutureValueSDP(ResSeg, WVTable, AlphaTable), stageprobResult)
  else
    #return  (FutureValueSDP(ResSeg, WVTable, AlphaTable, DValTable))
    return (FutureValueSDP(ResSeg, WVTable, AlphaTable), 0)
  end

end


# FUNCTIONS TO UPDATE STATE PARAMETERS BEFORE SOLVING STAGE PROBLEM
#---------------------------------------------------------------------------------------------

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
  ramping_constraints,
  #StepFlow,
 )
  Price = ScenarioLattice.states[t][iState, 2] .* PriceScale[t,1:NStep]                          # Per ogni stato (ce ne sono 5) della settimana corrente di quello scenario, moltiplico il prezzo per [0,75  1,0  1,25]
  
  for iMod = 1:HY.NMod                                                                           # Itero 2 volte (n. di bacini)  

    # set RHS for reservoir balance first time-step in stage 
    JuMP.set_normalized_rhs(
      SP.resbalInit[iMod],                                                                       # Per l'equazione "balance reservoir Initial state"
      ResSeg[nfrom][iMod] +                                                                      # Volume reservoir iMod nella cobinazione nfrom + fatt.conv* inflow(allo stato iState,settimana t)* scala 
      StepFranc[t,1] * ScenarioLattice.states[t][iState, 1] * HY.Scale[iMod],
    )
    #StepFranc[t,1]

    for iStep = 1:NStep                                                                          # per i 3 step
      # update power price to objective (per stage) 
      set_objective_coefficient(SP.model, SP.prod[iMod, iStep], NHoursStep * Price[iStep])       # Coefficiente che pongo davanti al prezzo nella funzione obiettivo
      set_objective_coefficient(SP.model, SP.pump[iMod, iStep], -NHoursStep * Price[iStep])  


      if iStep > 1
        # set RHS for reservoir balance time-step > 1
        JuMP.set_normalized_rhs(
          SP.resbalStep[iMod, iStep],
          StepFranc[t,iStep] * ScenarioLattice.states[t][iState, 1] * HY.Scale[iMod],
        )#StepFranc[t,iStep]
      end

      # MINIMUM ENVIRONMENTAL FLOW
      for n=1:(HY.N_min_flows[iMod]-1)                                                           # Cambio il valore di qMin: per determinate settimane ho valori diversi da 0
        if t>=HY.Activation_weeks[iMod,n] && t<HY.Activation_weeks[iMod,n+1]
            HY.qMin[iMod]= HY.Min_flows[iMod,n]
            JuMP.set_normalized_rhs(SP.q_min[iMod, iStep], HY.qMin[iMod])
        elseif t>=HY.Activation_weeks[iMod,n+1]
            HY.qMin[iMod]= HY.Min_flows[iMod,n+1]
            JuMP.set_normalized_rhs(SP.q_min[iMod, iStep], HY.qMin[iMod])
        end
      end

      # PUMP DEPEDENCE ON LOWER RESERVOIR LEVEL
      #if HY.Station_with_pump==1
       # SP = DeactivationPump_SDP(SP,t,HY,ResSeg,52,iStep,nfrom)
      #end

      # RAMPING CONSTRAINTS
      if iMod==1 && ramping_constraints
        SP = ramping_constraints_SDP(case::caseData,SP,iMod,t,nfrom,HY,ResSeg,iStep)
      end 

    end

  end
  return SP
end

function update_future_value_constraint(SP, AlphaTable, NSeg, t, iState)

  if HY.NMod == 2
    for iSegL = 1:NSeg[2]                                                    # iSegL=1:NSeg[2]
      for iSegU = 1:NSeg[1]                                                  # iSegU = 1:NSeg[1]
        JuMP.set_normalized_coefficient(
          SP.AlphaCon,
          SP.gamma[iSegU, iSegL],
          -AlphaTable[t+1, iState, (NSeg[1]*(iSegL-1)+iSegU)], #NSeg              #NSeg[1]*(iSegL-1)+iSegU
        )
      end
    end
  elseif HY.NMod == 1
    for iSeg = 1:NSeg[1]     #NSeg[1]
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

# FUNCTIONS TO UPDATE PARAMETERS OF END OF ITERATION
#---------------------------------------------------------------------------------------------

function check_diff(Table)
  if size(Table)[end] == 1                                                 # case of one reservoir - shouldn't put 1 instead of last :?
    diff_matrix = Table[1, :, :, :] .- Table[end, :, :, :]
    #diff = [sum(abs.(diff_matrix))]
    diff = abs.(diff_matrix)
  elseif size(Table)[end] == 2                                             # Case of 2 reservoirs - difference for upper and lower
    diff_u = Table[1, :, :, 1] .- Table[end, :, :, 1]
    diff_l = Table[1, :, :, 2] .- Table[end, :, :, 2]
    #diff = [sum(abs.(diff_u)), sum(abs.(diff_l))]
    diff = [abs.(diff_u); abs.(diff_l)]
  end
  return diff
end

function calculateWatervalues(nfrom, ResSeg, AlphaTable)      #nfrom = combo di volumi

  Table = zeros(HY.NMod)

  if HY.NMod == 2                                                                                                 # If we have 2 reservoirs
    if nfrom == 1                                                                                                 # Consider comb.n.1 = (0,00  0,00) - null volumes for both reservoirs
      Table[:] .= 0                                                                                               # WV in that point is nulla
    elseif nfrom <= NSeg[1]    #NSeg[1]     - first column -low reservoir is empty                                # If the upper reservoir gets filled but lower is empty
      Table[1] =                                                                                                  # WV of upper reservoir are updated
        (AlphaTable[nfrom] - AlphaTable[nfrom-1]) / (ResSeg[nfrom][1] - ResSeg[nfrom-1][1])
      Table[2] = 0                                                                                                # WV of lower reservoirs are null because empty
    elseif ((nfrom - 1) / NSeg[1]) % 1 == 0              #when upper reservoir is empty                           # In case when upper reservoi is empty and lower one is filling up
      Table[1] = 0                                                                                                # WV for upper reservoir are null    
      Table[2] =
        (AlphaTable[nfrom] - AlphaTable[nfrom-NSeg[1]]) /  #NSeg[1]
        (ResSeg[nfrom][2] - ResSeg[nfrom-NSeg[1]][2])        #NSeg[1]                                             # WV for lower reservoir are updated
    else
      Table[1] =
        (AlphaTable[nfrom] - AlphaTable[nfrom-1]) / (ResSeg[nfrom][1] - ResSeg[nfrom-1][1])                       # Otherwise , if both Upperres and Lower res have volumes >=0
      Table[2] =
        (AlphaTable[nfrom] - AlphaTable[nfrom-NSeg[1]]) /      #NSeg[1]
        (ResSeg[nfrom][2] - ResSeg[nfrom-NSeg[1]][2])          #NSeg[1]
    end
  elseif HY.NMod == 1                                                                                             # In the case we have only one reservoir  
    iMod = 1
    if nfrom == 1
      Table[iMod] = 0
    else
      Table[iMod] =
        (AlphaTable[nfrom] - AlphaTable[nfrom-1]) /
        (ResSeg[nfrom][iMod] - ResSeg[nfrom-1][iMod])
    end
  else
    error(" Number of modules not supported: ", HY.NMod)
  end
  return Table
end

function update_futureValue(WVTable, ResSeg, NEnvStates)

  Table = zeros(NStates * NEnvStates, length(ResSeg))       #update AlphaTable                                  # Creo una matrice nulla (5x25)

  if size(WVTable)[end] == 2                                                                                    # In case we have 2 reservoirs
    for iSeg = 2:length(ResSeg)
      segLower = Int(ceil(iSeg / NSeg[1])) #column number            #NSeg[2])
      segUpper = iSeg - (segLower - 1) * NSeg[2] #row number         #NSeg[1]
      if segLower == 1      #column 1 (lower reservoir is null)
        Table[:, iSeg] .=
          Table[:, iSeg-1] .+
          (ResSeg[iSeg][1] - ResSeg[iSeg-1][1]) .* WVTable[1, :, iSeg, 1]
      elseif segUpper == 1
        Table[:, iSeg] .=                                                         
          Table[:, iSeg-NSeg[1]] .+                                                      #NSeg[2] prima
          (ResSeg[iSeg][2] - ResSeg[iSeg-NSeg[1]][2]) .* WVTable[1, :, iSeg, 2]          #NSeg[1]
      else
        Table[:, iSeg] .=
          Table[:, iSeg-1] .+
          (ResSeg[iSeg][1] - ResSeg[iSeg-1][1]) .* WVTable[1, :, iSeg, 1]
      end
    end
  elseif size(WVTable)[end] == 1
    iMod = 1
    for iSeg = 2:length(ResSeg)
      Table[:, iSeg] .=
        Table[:, iSeg-1] .+
        (ResSeg[iSeg][iMod] - ResSeg[iSeg-1][iMod]) .* WVTable[1, :, iSeg, iMod]
    end
  end

  return Table
end


function isNotConvex(WVTable)
  notConvex = false
  if HY.NMod == 1
    for iSeg = 3:length(ResSeg)
      if any(WVTable[iSeg-1, :] .< WVTable[iSeg, :])
        notConvex = true
        break
      end
    end
  elseif HY.NMod == 2
    for iSeg = 2:length(ResSeg)
      if iSeg > NSeg[1]              #NSeg[1]

        if WVTable[iSeg-NSeg[1], 1] != 0                                       #NSeg[1]
          if (any(WVTable[iSeg-NSeg[1], 1] .< WVTable[iSeg, 1]))               #NSeg[1]
            notConvex = true
            break
          end
        end                                                                    # Lower Reservoir
        if WVTable[iSeg-NSeg[1], 2] != 0                                       # NSeg[1]
          if (any(WVTable[iSeg-NSeg[1], 2] .< WVTable[iSeg, 2]))               # NSeg[1]
            notConvex = true
            break
          end
        end
      elseif (((iSeg - 1) / NSeg[1]) % 1 != 0)                                  # NSeg[1]
        if WVTable[iSeg-1, 1] != 0
          if (any(WVTable[iSeg-1, 1] .< WVTable[iSeg, 1]))
            notConvex = true
            break
          end
        end
        if WVTable[iSeg-1, 2] != 0
          if (any(WVTable[iSeg-1, 2] .< WVTable[iSeg, 2]))
            notConvex = true
            break
          end
        end
      end
    end
  end
  return notConvex
end
