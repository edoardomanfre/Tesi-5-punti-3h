

function BuildProblem_sim(InputParameters::InputParam, HY::HydroData, SolverParameters::SolverParam)                  # A seconda che ho 1 o 2 reservoirs, vado alla relativa funzione per creare il MODELLO
  if HY.NMod == 1
    SP = BuildStageProblem_sim(InputParameters, HY, SolverParameters)
  elseif HY.NMod == 2
    SP = BuildStageProblemTwoRes_sim(InputParameters, HY, SolverParameters)
  end
  return SP
end

function BuildStageProblemTwoRes_sim(InputParameters::InputParam, HY::HydroData, SolverParameters::SolverParam)       # When we have 2 hydropower plants- 2 turbines

  @unpack (
    CPX_PARAM_SCRIND,
    CPX_PARAM_PREIND,
    CPXPARAM_MIP_Tolerances_MIPGap,
    CPX_PARAM_TILIM,
    CPX_PARAM_THREADS,
  ) = SolverParameters

  @unpack (NSeg, NStep, PenSpi, PenQ, AlphaMax, NStep, MM3Step, Big) = InputParameters

  M = Model(
    with_optimizer(
      CPLEX.Optimizer,
      CPX_PARAM_SCRIND = CPX_PARAM_SCRIND,
      CPX_PARAM_PREIND = CPX_PARAM_PREIND,
      CPXPARAM_MIP_Tolerances_MIPGap = CPXPARAM_MIP_Tolerances_MIPGap,
      CPX_PARAM_TILIM = CPX_PARAM_TILIM,
      CPX_PARAM_THREADS = CPX_PARAM_THREADS,
    ),
  )


  @variable(M,0 <= res[iMod = 1:HY.NMod, iStep = 1:NStep] <= HY.MaxRes[iMod], base_name = "res") #Inserire volume minimo diverso da 0

  @variable(M, 0 <= spi[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "spi")
  @variable(M, 0 <= prod[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "prod")

  @variable(M, 0 <= pump[iMod = 1:HY.NMod, iStep = 1:NStep]<= Big, base_name ="pump")                                         #-> pump is only present in the lower reservoir (n.2) - l'energia richiesta per il pompaggio non pu' essere negativa
                                                                                
  @variable(M, 0 <= q_slack[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "q_slack")
  @variable(M, 0 <= min_slack[iMod=1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "min_slack")

  @variable(M, 0 <= res_slack_pos[iMod=1:HY.NMod, iStep = 1:NStep] <= Big, base_name =" res_slack_pos")
  @variable(M, 0 <= res_slack_neg[iMod=1:HY.NMod, iStep = 1:NStep] <= Big, base_name =" res_slack_neg")

  # variabile per caratteristiche della turbina
  @variable(M,0 <=disSeg[iMod = 1:HY.NMod, iSeg = 1: (HY.NDSeg[iMod]-1), iStep = 1:NStep] <= HY.DisPointTurb[iMod, iSeg+1],base_name = "dseg")

  #Variable for power turbine 
  @variable(M,0<=powSegTurb[iMod = 1:HY.NMod, iSeg = 1: (HY.NDSeg[iMod]-1), iStep = 1:NStep]<= HY.PowMaxSegTurb[iMod, iSeg+1],base_name ="psegturb")

  #Binary variable for turbine
  @variable(M,u_turb_1[iMod = 1:HY.NMod, iStep = 1:NStep], Bin, base_name = "bin_turb_1")
  @variable(M,u_turb_2[iMod = 1:HY.NMod, iStep = 1:NStep], Bin, base_name = "bin_turb_2")
  @variable(M,u_turb_3[iMod = 1:HY.NMod, iStep = 1:NStep], Bin, base_name = "bin_turb_3")
  @variable(M,u_turb_4[iMod = 1:HY.NMod, iStep = 1:NStep], Bin, base_name = "bin_turb_4")

  # Variable for pumping discharge
  @variable(M,0<=disSegPump[iStep = 1:NStep]<= HY.DisPointPump[HY.NDSegPump],base_name ="dsegpump")

  #Variable for power pump
  @variable(M,0<=powSegPump[iStep = 1:NStep]<= HY.PowMaxSegPump[HY.NDSegPump],base_name ="psegpump")

  #Binary variable for pump
  @variable(M,u_pump[iStep = 1:NStep], Bin, base_name = "bin_pump")

  # Variabile By-pass per deflusso minimo vitale
  @variable(M,by_pass[iMod = 1:HY.NMod,iStep=1:NStep]>=0, base_name = "min_flow_bypass")


  @variable(M, alpha <= Big, base_name = "alp")
  @variable(M, 0 <= beta_upper[nU = 1:NSeg[1]] <= 1, base_name = "weightReservoirUpper")
  @variable(M, 0 <= beta_lower[nL = 1:NSeg[2]] <= 1, base_name ="weightReservoirLower")
  #@variable(M, 0 <= χ[i = 1:(2*NSeg-1)] <= 1, base_name = "weightReservoirDiag")
  @variable(M, 0 <= gamma[nU = 1:NSeg[1], nL = 1:NSeg[2]] <= 1, base_name = "weight")       #nL = 1:NSeg*7


 #CHANGE OBJECTIVE function
 @objective(
  M,
  MathOptInterface.MAX_SENSE,
  sum(
    1*prod[iMod,iStep]-1*pump[iMod,iStep]- 
    PenSpi * spi[iMod,iStep] - 
    1E5 * q_slack[iMod,iStep] - 1E5 * min_slack[iMod,iStep] - 
    1E5*res_slack_pos[iMod,iStep] -1E5*res_slack_neg[iMod,iStep] for iMod = 1:HY.NMod for iStep = 1:NStep)+
    alpha ) 



  #Reservoir balances for first step considering PUMPING 
  @constraint(
    M,resbalInit[iMod = 1:HY.NMod],         
    res[iMod, 1] +                                                                                                                                                                                                                                             
    MM3Step * sum(disSeg[iMod, iSeg, 1] for iSeg = 1:(HY.NDSeg[iMod]-1)) +     
    MM3Step * spi[iMod, 1] -                                                                                                
    MM3Step * sum(disSeg[uMod, iSeg, 1] for uMod=1:HY.NUp[iMod] for iSeg = 1:(HY.NDSeg[uMod]-1)) -
    MM3Step * sum(spi[uMod,1] for uMod=1:HY.NUp[iMod])-                                                                                                
    MM3Step * disSegPump[1] * HY.Pump_direction[iMod] + 
    MM3Step * by_pass[iMod,1] -
    MM3Step * sum(by_pass[uMod,1] for uMod=1:HY.NUp[iMod]) == 0                                                                                                                                       
  ) #resbalInitRHS[iMod] + StepFranc*inf[iMod])

#Reservoir balance within stage with PUMPING
@constraint(
  M,resbalStep[iMod=1:HY.NMod, iStep = 2:NStep],
  res[iMod, iStep] - res[iMod, iStep-1] +
  MM3Step * sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:(HY.NDSeg[iMod]-1)) +
  MM3Step * spi[iMod, iStep] -
  MM3Step * sum(disSeg[uMod, iSeg, iStep] for uMod=1:HY.NUp[iMod] for iSeg = 1:(HY.NDSeg[uMod]-1)) -
  MM3Step * sum(spi[uMod,iStep] for uMod=1:HY.NUp[iMod])- 
  MM3Step * (disSegPump[iStep]*HY.Pump_direction[iMod]) +
  MM3Step * by_pass[iMod,iStep] -
  MM3Step * sum(by_pass[uMod,iStep] for uMod=1:HY.NUp[iMod])  == 0
)


  # Lower reservoir volume constraint (slack with punishment)
  @constraint(
    M,
    minResPunish[iMod = 1:HY.NMod, iStep = 1:NStep],
    res[iMod, iStep] >= HY.MaxRes[iMod] * 0
  )


  #Hydropower generation
  @constraint(
    M,
    prodeff[iMod = 1:HY.NMod, iStep = 1:NStep],
    prod[iMod, iStep] == sum(powSegTurb[iMod, iSeg, iStep] for iSeg = 1:(HY.NDSeg[iMod]-1))
  )
  
  #=
@constraint(
    M,
    pumpdischarge[tSeg=1:HY.NDSegPump,iStep=1:NStep],
    disSegPump[tSeg, iStep] <= HY.DisMaxSegPump[tSeg]
  ) =#

  #Pumping power required
  @constraint(
    M,
    pumpeff_up[iMod=1, iStep=1:NStep],
    pump[iMod,iStep] == powSegPump[iStep]
  )

  @constraint(
    M,
    pumpeff_low[iMod=2, iStep=1:NStep],
    pump[iMod,iStep] == 0
  )

  # Minimum environmental flow
  @constraint(
    M,
    q_min[iMod = 1:HY.NMod, iStep=1:NStep], 
    by_pass[iMod,iStep] + min_slack[iMod,iStep] == HY.qMin[iMod])   #>

  # Ramping constraints
 @constraint(
    M,
    positive_var[iMod=1:HY.NMod,iStep=2:NStep],
    res[iMod,iStep] - res_slack_pos[iMod,iStep]<= res[iMod,iStep-1] + 1000   # -res_slack_pos[iMod,iStep] 
    )
  
  @constraint(
    M,
    negative_var[iMod=1:HY.NMod,iStep=2:NStep],
    res[iMod,iStep] + res_slack_neg[iMod,iStep]>= res[iMod,iStep-1] - 1000  # +res_slack_neg[iMod,iStep]
    )

 @constraint(
      M,
      Initial_volumevar_positive[iMod=1:HY.NMod],
      res[iMod,1] - res_slack_pos[iMod,1]  <= HY.MaxRes[iMod] 
      )
  
  @constraint(
      M,
      Initial_volumevar_negative[iMod=1:HY.NMod],
      res[iMod,1] + res_slack_neg[iMod,1]  >= 0 
      )

#SOS2 constraint
  @constraint(
    M,
    AlphaCon,
    alpha - sum(gamma[iSegU, iSegL] for iSegU = 1:NSeg[1] for iSegL = 1:NSeg[2]) == 0         #iSegL=1:NSeg*7
  ) #Also updated in SDP

  @constraint(
    M,
    weightcon_upper[iMod = 1, nU = 1:NSeg[1]],
    beta_upper[nU] == sum(gamma[nU, nL] for nL = 1:NSeg[2])                                # nL = 1: NSeg*7
  )
  @constraint(
    M,
    weightcon_lower[iMod = 2, nL = 1:NSeg[2]],                                             # nL = 1:NSeg*7
    beta_lower[nL] == sum(gamma[nU, nL] for nU = 1:NSeg[1])
  )

  #End reservoir in last time step, using SOS2 weights
  @constraint(
    M,
    endResWeight_upper[iMod = 1],
    res[iMod, end] == sum(beta_upper[iSeg] * ResSeg[iSeg, 1][1] for iSeg = 1:NSeg[1])
  )
  @constraint(
    M,
    endResWeight_lower[iMod = 2],
    res[iMod, end] == sum(beta_lower[iSeg] * ResSeg[1, iSeg][2] for iSeg = 1:NSeg[2])          # iSeg=1:NSeg*7
  )

  #Binding the weights to SOS2 and the sum of the weights to be 1
  @constraint(
    M,
    weightcon,
    sum(gamma[iSegU, iSegL] for iSegU = 1:NSeg[1] for iSegL = 1:NSeg[2]) == 1           # iSeg=1:NSeg*7
  )
  #@constraint(M, SOS2[iMod=1:HY.NMod], beta[iMod,:] in MOI.SOS2(collect(1:NSeg)))

  @constraint(
    M,
    maxRelease[iMod = 1:HY.NMod, iStep = 1:NStep],
    sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:(HY.NDSeg[iMod]-1)) <= HY.DisPointTurb[iMod, HY.NDSeg[iMod]]
  )

 @constraint(
    M,
    maxReleasePump[iStep = 1:NStep],
    disSegPump[iStep] <= HY.DisPointPump[HY.NDSegPump]
  )

  #Pump first segment 
  @constraint(
    M,
    l_p_seg[iStep = 1:NStep],
    disSegPump[iStep] >= HY.DisPointPump[1]*u_pump[iStep]
  )

  @constraint(
    M,
    u_p_seg[iStep = 1:NStep],
    disSegPump[iStep] <= HY.DisPointPump[2]*u_pump[iStep]
  )

  @constraint(
    M,
    maxPowerPump[iStep = 1:NStep],
    powSegPump[iStep] == disSegPump[iStep] * HY.m_pump + u_pump[iStep] * HY.K_pump_max
  )

  #Turbine first segment
  @constraint(
    M,
    l_t_seg_1[iMod = 1:HY.NMod, iSeg = 1, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] >= HY.DisPointTurb[iMod, 1]*u_turb_1[iMod, iStep]
  )

  @constraint(
    M,
    u_t_seg_1[iMod = 1:HY.NMod, iSeg = 1, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] <= HY.DisPointTurb[iMod, 2]*u_turb_1[iMod, iStep]
  )

  @constraint(
    M,
    maxPowerTurb_1[iMod = 1:HY.NMod, iSeg = 1, iStep = 1:NStep],
    powSegTurb[iMod, iSeg, iStep] == disSeg[iMod, iSeg, iStep] * HY.m_1[iMod]  + u_turb_1[iMod, iStep] * HY.K_1_max[iMod]
  )

  #Turbine second segment
  @constraint(
    M,
    l_t_seg_2[iMod = 1:HY.NMod, iSeg = 2, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] >= HY.DisPointTurb[iMod, 2]*u_turb_2[iMod, iStep]
  )

  @constraint(
    M,
    u_t_seg_2[iMod = 1:HY.NMod, iSeg = 2, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] <= HY.DisPointTurb[iMod, 3]*u_turb_2[iMod, iStep]
  )

  @constraint(
    M,
    maxPowerTurb_2[iMod = 1:HY.NMod, iSeg = 2, iStep = 1:NStep],
    powSegTurb[iMod, iSeg, iStep] == disSeg[iMod, iSeg, iStep] * HY.m_2[iMod] + u_turb_2[iMod, iStep] * HY.K_2_max[iMod]
  )

  #Turbine third segment
  @constraint(
    M,
    l_t_seg_3[iMod = 1:HY.NMod, iSeg = 3, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] >= HY.DisPointTurb[iMod, 3]*u_turb_3[iMod, iStep]
  )

  @constraint(
    M,
    u_t_seg_3[iMod = 1:HY.NMod, iSeg = 3, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] <= HY.DisPointTurb[iMod, 4]*u_turb_3[iMod, iStep]
  )

  @constraint(
    M,
    maxPowerTurb_3[iMod = 1:HY.NMod, iSeg = 3, iStep = 1:NStep],
    powSegTurb[iMod, iSeg, iStep] == disSeg[iMod, iSeg, iStep] * HY.m_3[iMod] + u_turb_3[iMod, iStep] * HY.K_3_max[iMod]
  )

  #Turbine four segment
  @constraint(
    M,
    l_t_seg_4[iMod = 1:HY.NMod, iSeg = 4, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] >= HY.DisPointTurb[iMod, 4]*u_turb_4[iMod, iStep]
  )

  @constraint(
    M,
    u_t_seg_4[iMod = 1:HY.NMod, iSeg = 4, iStep = 1:NStep],
    disSeg[iMod, iSeg, iStep] <= HY.DisPointTurb[iMod, 5]*u_turb_4[iMod, iStep]
  )

  @constraint(
    M,
    maxPowerTurb_4[iMod = 1:HY.NMod, iSeg = 4, iStep = 1:NStep],
    powSegTurb[iMod, iSeg, iStep] == disSeg[iMod, iSeg, iStep] * HY.m_4[iMod] + u_turb_4[iMod, iStep] * HY.K_4_max[iMod]
  )

  @constraint(
    M,
    binaries[iMod = 1:HY.NMod, iStep = 1:NStep],
    u_turb_1[iMod, iStep] + u_turb_2[iMod, iStep] + u_turb_3[iMod, iStep] + u_turb_4[iMod, iStep] + u_pump[iStep] <= 1
  )

  @constraint(M, minReservoirEnd[iMod = 1:HY.NMod, iStep = NStep], res[iMod, iStep] >= 0)

  @constraint(M, minReservoir[iMod = 1:HY.NMod, iStep = 1:NStep], res[iMod, iStep] >= 0)

  @constraint(M, noDecrease_week[iMod = 1:HY.NMod, iStep = NStep], res[iMod, iStep] >= 0)

  return StageProblemTwoRes_sim(
    M,
    res,
    spi,
    prod,
    pump,
    q_slack,
    min_slack,
    res_slack_pos,
    res_slack_neg,
    q_min,
    disSeg,
    disSegPump,
    powSegPump,
    u_pump,
    powSegTurb,
    u_turb_1,
    u_turb_2,
    u_turb_3,
    u_turb_4,
    by_pass,
    resbalInit,
    resbalStep,
    positive_var,
    negative_var,
    Initial_volumevar_positive,
    Initial_volumevar_negative,
    prodeff,
    #pumpdischarge,
    pumpeff_up,
    pumpeff_low,
    alpha,
    gamma,
    AlphaCon,
    beta_upper,
    beta_lower,
    #χ,
    maxRelease,
    maxReleasePump,
    maxPowerPump,
    l_p_seg,
    u_p_seg,
    maxPowerTurb_1,
    l_t_seg_1,
    u_t_seg_1,
    maxPowerTurb_2,
    l_t_seg_2,
    u_t_seg_2,
    maxPowerTurb_3,
    l_t_seg_3,
    u_t_seg_3,
    maxPowerTurb_4,
    l_t_seg_4,
    u_t_seg_4,
    binaries,
    minReservoirEnd,
    minReservoir,
    noDecrease_week,
    minResPunish,
  )
end
    #q_min,

function BuildStageProblem_sim(InputParameters::InputParam, HY::HydroData, SolverParameters::SolverParam)                                     # Used only for 1 hydropower plant - 1 turbine
    @unpack (
      CPX_PARAM_SCRIND,
      CPX_PARAM_PREIND,
      CPXPARAM_MIP_Tolerances_MIPGap,
      CPX_PARAM_TILIM,
      CPX_PARAM_THREADS,
    ) = SolverParameters
  
    @unpack (NSeg, NStep, PenSpi, PenQ, AlphaMax, NStep, MM3Step, Big) = InputParameters

  M = Model(
    with_optimizer(
      CPLEX.Optimizer,
      CPX_PARAM_SCRIND = CPX_PARAM_SCRIND,
      CPX_PARAM_PREIND = CPX_PARAM_PREIND,
      CPXPARAM_MIP_Tolerances_MIPGap = CPXPARAM_MIP_Tolerances_MIPGap,
      CPX_PARAM_THREADS = 1
    ),
  )

  @variable( M, 0 <= res[iMod = 1:HY.NMod, iStep = 1:NStep] <= HY.MaxRes[iMod],base_name = "res")

  @variable(M, 0 <= spi[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "spi")
  @variable(M, 0 <= prod[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "prod")
  @variable(M, 0 <= q_slack[iMod = 1:HY.NMod, iStep = 1:NStep] <= Big, base_name = "q_slack")
  @variable(M, 0 <= min_slack[iMod=1:HY.NMod, iStep=1:NStep] <= Big, base_name="min_slack")

  @variable(M,0 <= disSeg[iMod = 1:HY.NMod, iSeg = 1:(HY.NDSeg[iMod]-1), iStep = 1:NStep] <= HY.DisMaxSeg[iMod, iSeg],base_name = "dseg")
  @variable(M,by_pass[iMod = 1:HY.NMod,iStep=1:NStep]>=0, base_name = "min_flow_bypass")

  @variable(M, alpha <= Big, base_name = "alp")
  @variable(M, 0 <= gamma[n = 1:NSeg[1]] <= 1, base_name = "weight")

  @objective(
    M,
    MathOptInterface.MAX_SENSE,
    sum(
      1*prod[iMod, iStep] - PenSpi * spi[iMod, iStep] - q_slack[iMod, iStep] * 10E5 - 1E5 * min_slack[iMod,iStep] for iMod = 1:HY.NMod for iStep = 1:NStep) + alpha
  )

  #Reservoir balances for first step
  @constraint(
    M,
    resbalInit[iMod = 1:HY.NMod],
    res[iMod, 1] +
    #MM3Step * sum(disSeg[iMod, iSeg, 1] for iSeg = 1:HY.NDSeg[iMod]) +
    MM3Step * sum(disSeg[iMod, iSeg, 1] for iSeg = 1:(HY.NDSeg[iMod]-1)) +
    MM3Step * spi[iMod, 1] -
    #MM3Step * sum(disSeg[uMod, iSeg, 1] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NDSeg[uMod]) -
    MM3Step * sum(disSeg[uMod, iSeg, 1] for uMod = 1:HY.NUp[iMod] for iSeg = 1:(HY.NDSeg[uMod]-1)) -
    MM3Step * sum(spi[uMod, 1] for uMod = 1:HY.NUp[iMod]) +
    MM3Step * by_pass[iMod,1] -
    MM3Step * sum(by_pass[uMod,1] for uMod=1:HY.NUp[iMod])  == 0 
  )#resbalInitRHS[iMod] + StepFranc*inf[iMod])

  #Reservoir balance within stage
  @constraint(
    M,
    resbalStep[iMod = 1:HY.NMod, iStep = 2:NStep],
    res[iMod, iStep] - res[iMod, iStep-1] +
#    MM3Step * sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:HY.NDSeg[iMod]) +
    MM3Step * sum(disSeg[iMod, iSeg, iStep] for iSeg = 1:(HY.NDSeg[iMod]-1)) +
    MM3Step * spi[iMod, iStep] -
#    MM3Step * sum(disSeg[uMod, iSeg, iStep] for uMod = 1:HY.NUp[iMod] for iSeg = 1:HY.NDSeg[uMod]) -
    MM3Step * sum(disSeg[uMod, iSeg, iStep] for uMod = 1:HY.NUp[iMod] for iSeg = 1:(HY.NDSeg[uMod]-1)) -
    MM3Step * sum(spi[uMod, iStep] for uMod = 1:HY.NUp[iMod]) +
    MM3Step * by_pass[iMod,iStep] -
    MM3Step * sum(by_pass[uMod,iStep] for uMod=1:HY.NUp[iMod]) == 0
  ) #StepFranc*inf[iMod])

  # Lower reservoir volume constraint (slack with punishment)
  @constraint(
    M,
    minResPunish[iMod = 1:HY.NMod, iStep = 1:NStep],
    res[iMod, iStep] >= HY.MaxRes[iMod] * 0
  )#0.1)

  #Hydropower generation
  @constraint(
    M,
    prodeff[iMod = 1:HY.NMod, iStep = 1:NStep],
    prod[iMod, iStep] == powSegTurb[iMod, iSeg, iStep]
  )

  # Min flow requirment
  @constraint(M,q_min[iMod = 1:HY.NMod, iStep=1:NStep], by_pass[iMod,iStep] + min_slack[iMod,iStep] >= HY.qMin[iMod])

  # RAMPING CONSTRAINTS
  @constraint(
    M,
    positive_var[iMod=1:HY.NMod,iStep=2:NStep],
    res[iMod,iStep]<= res[iMod,iStep-1] + 1000
    )
  
  @constraint(
    M,
    negative_var[iMod=1:HY.NMod,iStep=2:NStep],
    res[iMod,iStep]>= res[iMod,iStep-1] -1000
    )

 @constraint(
      M,
      Initial_volumevar_positive[iMod=1:HY.NMod],
      res[iMod,1]<= HY.MaxRes[iMod]
      )
  
  @constraint(
      M,
      Initial_volumevar_negative[iMod=1:HY.NMod],
      res[iMod,1]>= 0
      )

  @constraint(M, AlphaCon, alpha - sum(gamma[iSeg] for iSeg = 1:NSeg[1]) == 0)

  @constraint(
    M,
    endResWeight[iMod = 1:HY.NMod],
    res[end] == sum(gamma[iSeg] * ResSeg[iSeg][iMod] for iSeg = 1:NSeg[1])
  )

  @constraint(M, 
    weightcon, 
    sum(gamma[iSeg] for iSeg = 1:NSeg[1]) == 1)

  @constraint(
    M,
    maxRelease[iMod = 1:HY.NMod, iStep = 1:NStep],
    disSeg[iMod, iStep] <= HY.DisMaxSeg[iMod]
  )

  @constraint(M, minReservoirEnd[iMod = 1:HY.NMod, iStep = NStep], res[iMod, iStep] >= 0)

  @constraint(M, minReservoir[iMod = 1:HY.NMod, iStep = 1:NStep], res[iMod, iStep] >= 0)

  @constraint(M, noDecrease_week[iMod = 1:HY.NMod, iStep = NStep], res[iMod, iStep] >= 0)

  return StageProblem_sim(
    M,
    res,
    spi,
    prod,
    q_slack,
    min_slack,
    q_min,
    disSeg,
    by_pass,
    resbalInit,
    resbalStep,
    positive_var,
    negative_var,
    Initial_volumevar_positive,
    Initial_volumevar_negative,
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

