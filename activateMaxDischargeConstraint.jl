

# check and activate environmental constraint, update problem                             # ma e' solo per il volume Lower su cui abbaimo i vincoli
#------------------------------------------------------------------

function activate_EnvConstraint(                                                          # Chiede in input : il modello, la settimana t, lo stato del sistema in tale settimana, la cobinazione dei due volumi , il reticolo degli scenari,
  SP,                                                                                     # i dati idraulici dei bacini , deflusso minimo vitale, n.step nella settimana, envData = boolean (false/true) se il limite e' stato passato o no
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
  pump_cons = add_dischargeLimitPump
  if t >= envData.lastAct && t <= envData.lastMaxDisch                    #se siamo tra settimana 23 e 38: se t>= last activation week(23) e t<= last_active_wek_max_discharge(38)                                                               
    SP = relax_noResDecrease(SP, envData.envMod, NStep)                   #Volume[2,NStep_ultimo]>=0, può diminuire quanto vuole
    
    if ResSeg[nfrom][envData.envMod] >= envData.deactLevel                #se in quelle settimane (23-38), il volume del bacino su cui è posto il vincolo >= livello di deactivation (87.44)
      SP = add_minRes(SP, envData, NStep)                                 #allora il deactLevel (87.44 Mmn3) diventa il volume minimo per tutta la settimana
#      add_dischargeLimitPump = false                                      #non ho limiti sul pompaggio se non vado sotto il deactLevel
      pump_cons = false
      #SP = relax_minRes_init(SP,envData)
      #SP = relax_dichargeLimit(SP,envData)
    elseif (                                                              #se volume del bacino Imod allo stato nfrom + (inflow allo stato iState della settimana t)*fattore_scala(bacino iMod) >= 87.44 Mm^3
      (ResSeg[nfrom][envData.envMod] + ScenarioLattice.states[t][iState, 1] * HY.Scale[envData.envMod]) >= envData.deactLevel
    )
      SP = add_minRes_init(SP, envData, NStep)                            #se il volume del bacino+incoming inflow >= deactLevel -> allora deactLevel diventa volume minimo                                       
#      add_dischargeLimitPump = false                                      #non ho limiti sul pompaggio
      pump_cons = false
      #SP = relax_dichargeLimit(SP,envData)
      #SP = relax_minRes(SP,envData)
    else                                                                                                                        
      SP = add_dichargeLimit(SP, envData, NStep)   
#      add_dischargeLimitPump = true                                       #ho limiti sul pompaggio: non posso scendere sotto 52 Mm3             
      pump_cons = true                    
      #SP = add_dischargeLimitPump(SP,NStep)
      #SP = relax_minRes_init(SP,envData)
      #SP = relax_minRes(SP,envData)
    end
  
  elseif t >= envData.lastAct && t <= envData.lastNoDecrease                              # Se t>=23 e t<=38
    SP = add_noResDecrease(SP, envData, ResSeg[nfrom][envData.envMod], NStep)             # non posso diminuire il volume d'acqua                 
  end

  return SP, pump_cons
end

function activate_EnvConstraint_sim(                                                      #funzione applicata nel codice "simulazione" con 100 scenari scelti
  SP,
  t,
  iScen,
  scenarios,
  HY,
  Reservoir,
  earlyActive_maxDischarge,
  envData,
  NStep,
  add_dischargeLimitPump,  #inizia con false
)
  pump_cons = add_dischargeLimitPump
  if envData.firstAct != 0 && t >= envData.firstAct && t < envData.lastAct                #not present 0 first - 23 last
    
    if earlyActive_maxDischarge                                                       
      SP = add_dichargeLimit(SP, envData, NStep)
#      add_dischargeLimitPump = true
      pump_cons = true
    elseif scenarios[iScen][t, 1] * HY.Scale[envData.envMod] >= envData.actLevel
      earlyActive_maxDischarge = true
    end
    
  elseif t >= envData.lastAct && t <= envData.lastMaxDisch        # for weeks between 23 and 38       87.44
    
    if Reservoir[envData.envMod, iScen, t-1, end] >= envData.deactLevel 
      SP = add_minRes(SP, envData, NStep)   # aggiungo livello minimo del Volume (non può andare sotto)
      SP = relax_minRes_init(SP, envData.envMod, NStep) # a fine settimana, posso andare anche sotto
      SP = relax_dichargeLimit(SP, envData.envMod, NStep) #scarico as much as I can
#      add_dischargeLimitPump = false
      pump_cons = false
    elseif (Reservoir[envData.envMod, iScen, t-1, end] + scenarios[iScen][t, 1] * HY.Scale[envData.envMod]) >= envData.deactLevel
      SP = add_minRes_init(SP, envData, NStep)    #a fine settimana devo raggiungere il limite iStep=NStep: volume >= 87.44
      SP = relax_dichargeLimit(SP, envData.envMod, NStep) #posso turbinare
#      add_dischargeLimitPump = false        #non metto limiti sullo scarico della pompa
      pump_cons = false
    else
      SP = add_dichargeLimit(SP, envData,NStep)
#      add_dischargeLimitPump = true     #non permetto che scarichi dalla pompa
      pump_cons = true
    end

  elseif t > envData.lastMaxDisch && t <= envData.lastNoDecrease      # No weeks - non usato nel nostro caso
    SP = add_noResDecrease(SP, envData, Reservoir[envData.envMod, iScen, t-1, end], NStep)
    SP = relax_dichargeLimit(SP, envData.envMod, NStep)
#    add_dischargeLimitPump = false
    pump_cons = false
    SP = relax_minRes_init(SP, envData.envMod, NStep)
    SP = relax_minRes(SP, envData.envMod, NStep)

  elseif t > envData.lastNoDecrease                             # for t>38
    SP = relax_noResDecrease(SP, envData.envMod, NStep)
    SP = relax_dichargeLimit(SP, envData.envMod, NStep)
    SP = relax_minRes_init(SP, envData.envMod, NStep)
    SP = relax_minRes(SP, envData.envMod, NStep)
#    add_dischargeLimitPump = false      # non ho limiti sul pompaggio
    pump_cons = false
  end

  return SP, earlyActive_maxDischarge, pump_cons
end

# functions to update stageproblem in JuMP
#----------------------------------------------

# discharge limit
function add_dichargeLimit(SP, envData, NStep)                         # attivo i vincoli sullo scarico
  iMod = envData.envMod                                                # iMod=2 perche' il vincolo e' sul bacino inferiore
  for iStep = 1:NStep
    set_normalized_rhs(SP.maxRelease[envData.envMod, iStep], envData.maxDischarge)
  end
  return SP
end

function add_disLimitPump(SP, NStep)                                   # attivo i vincoli sullo scarico                                                                         
  for iStep = 1:NStep
    set_normalized_rhs(SP.maxReleasePump[iStep], 0)
  end
  return SP
end

#=
function add_dischargeLimitPump(SP,NStep)                              # attivo i vincoli sullo scarico                                                                           # iMod=2 perche' il vincolo e' sul bacino inferiore
  for iStep = 1:NStep
    set_normalized_rhs(SP.pumpdischarge[iStep],0)
  end
  return SP
end
=#

function relax_dichargeLimit(SP, iMod, NStep)                          # disattivo i vincoli sullo scarico
  #    iMod = envData.envMod
  for iStep = 1:NStep
    set_normalized_rhs(
      SP.maxRelease[iMod, iStep],
      sum(HY.DisMaxSeg[iMod, iSeg] for iSeg = 1:HY.NDSeg[iMod]),
    )
  end
  return SP
end

function relax_disLimitPump(SP,NStep)                                  # disattivo i vincoli sullo scarico
  #    iMod = envData.envMod
  for iStep = 1:NStep
    set_normalized_rhs(
      SP.maxReleasePump[iStep],
      sum(HY.DisMaxSegPump[tSeg] for tSeg = 1:HY.NDSegPump),
    )
  end
  return SP
end

#=
function relax_dischargeLimitPump(SP,NStep)                            # disattivo i vincoli sullo scarico
  for iStep = 1:NStep
    for tSeg=1:HY.NDSegPump
    set_normalized_rhs(
      SP.pumpdischarge[tSeg,iStep], 
      HY.DisMaxSegPump[tSeg] 
      )
    end
  end
  return SP
end
=#

# min. reservoir level
function add_minRes_init(SP, envData, NStep)
  iMod = envData.envMod
  iStep = NStep
  set_normalized_rhs(SP.minReservoirEnd[iMod, iStep], envData.deactLevel) #87.44
  return SP
end

function relax_minRes_init(SP, iMod, NStep)
  #  iMod = envData.envMod
  iStep = NStep
  set_normalized_rhs(SP.minReservoirEnd[iMod, iStep], 0)
  return SP
end

function add_minRes(SP, envData, NStep)
  iMod = envData.envMod
  for iStep = 1:NStep
    set_normalized_rhs(SP.minReservoir[iMod, iStep], envData.deactLevel)
  end
  return SP
end

function relax_minRes(SP, iMod, NStep)
  #  iMod = envData.envMod
  for iStep = 1:NStep
    set_normalized_rhs(SP.minReservoir[iMod, iStep], 0)
  end
  return SP
end

# decrease in reservoir level not permitted
function add_noResDecrease(SP, envData, resPrev, NStep)
  iMod = envData.envMod
  iStep = NStep
  set_normalized_rhs(SP.noDecrease_week[iMod, iStep], resPrev)
  #@constraint(SP.model, noDecrease_init[iMod, iStep=1],  SP.res[iMod,iStep] >=  ResSeg[nfrom][2])
  #@constraint(SP.model, noDecrease[iMod, iStep=2:NStep],  SP.res[iMod,iStep] >= SP.res[iMod,iStep-1])
  return SP
end

function relax_noResDecrease(SP, iMod, NStep)
  # iMod = envData.envMod
  iStep = NStep
  set_normalized_rhs(SP.noDecrease_week[iMod, iStep], 0)
  #@constraint(SP.model, noDecrease_init[iMod, iStep=1],  SP.res[iMod,iStep] >=  ResSeg[nfrom][2])
  #@constraint(SP.model, noDecrease[iMod, iStep=2:NStep],  SP.res[iMod,iStep] >= SP.res[iMod,iStep-1])
  return SP
end
