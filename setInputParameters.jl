



# General read and formatting functions
#--------------------------------------

function read_from_file_to_dict!(f, vars::Dict{Symbol,Bool})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = parse(Bool, val)
      catch
      end

    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,Float64})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = parse(Float64, val)
      catch
      end

    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,String})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = val
      catch
      end
    end
  end
  return vars
end

function read_from_file_to_dict!(f, vars::Dict{Symbol,Any})
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = val
      catch
      end
    end
  end
  return vars
end


function read_from_file_to_dict!(f, vars::Dict{Symbol,Number})
  # Use Dict{Symbol, Number} to allow float and int values
  # NB: all vaues read as float from file initially
  for line in eachline(f)
    if !isempty(strip(line)) && !startswith(strip(line), "#")
      var, val = strip.(split(line, "="))
      try
        vars[Symbol(strip(var))] = parse(Float64, val)
      catch
      end
    end
  end
  return vars
end

function read_type_to_dict(file, type::DataType)
  f = open(file)
  paramDict = Dict{Symbol,type}()
  paramDict = read_from_file_to_dict!(f, paramDict)

  return paramDict
end

function set_integers!(paramDict, integers)
  for (key, value) in paramDict
    if key in integers
      paramDict[key] = parse(Int64,value)
    end
  end
  return paramDict
end

function set_floats!(paramDict, float)
  for (key, value) in paramDict
    if key in float
      paramDict[key] = parse(Float64,value)
    end
  end
  return paramDict
end

function set_vector!(paramDict, vector)
  value = paramDict[vector]
  items = split(value," ")
  NSeg = zeros(Int64,2)
  NSeg[1] = parse(Int64,items[1])
  NSeg[2] = parse(Int64,items[2])
  return NSeg
end

function read_csv(file, path)
  try
    println("Reading from: ", file)
    data = Matrix{Float64}(CSV.read(joinpath(path, file), DataFrame, header = false)) #Mm3 and €/MWh
    return data
  catch e
    println("Can't read ", file)
      throw(error())
  end
end

# Functions to read and prep from predefined files
#-----------------------------------------

# parameters
function read_parameters_from_config_file(file = "configParameters.in")

  #paramDict = read_type_to_dict(file,Number)
  paramDict = read_type_to_dict(file, Any)
  println("Parameters to be used:", paramDict)

  integers =[:NPrice :NStage :NHoursStage :NStep :NStates :MaxIt :NSamples :NSimScen :CPX_PARAM_SCRIND :CPX_PARAM_PREIND :CPX_PARAM_TILIM :CPX_PARAM_THREADS :LimitPump]
  paramDict = set_integers!(paramDict, integers)

  floats =[:PenSpi :PenQ :AlphaMax :conv :Big]
  paramDict = set_floats!(paramDict,floats)

  vector = :NSeg
  NSeg = set_vector!(paramDict,vector)

  paramDict[:MM3Week] = 60 * 60 * paramDict[:NHoursStage] / 10E5 #conversion factor from m3/s to mm3 per stage,
  paramDict[:NHoursStep] = Int(paramDict[:NHoursStage] / paramDict[:NStep])
  paramDict[:MM3Step] = paramDict[:MM3Week] / paramDict[:NStep]
  #paramDict[:StepFlow] = paramDict[:NHoursStep] / paramDict[:NHoursStage]

  StepFranc = read_csv("Inflow_ita_3h.csv",case.DataPath)

  #integers =
   # [:NPrice :NStage :NSeg :NHoursStage :NHoursStep :NStep :NStates :MaxIt :NSamples :NSimScen :CPX_PARAM_SCRIND :CPX_PARAM_PREIND :CPX_PARAM_TILIM]
  #paramDict = set_integers!(paramDict, integers)

  #inputData = InputParam(; paramDict...,StepFranc)
  inputData = InputParam(;paramDict..., NSeg, StepFranc)

  return inputData
end

function read_solverParameters_from_file(file = "solverParameters.in")
  try
    println("Reading solver parameters from config file.")
    #paramDict = read_type_to_dict(file, Number)
    paramDict = read_type_to_dict(file, Any)

    #integers = [:CPX_PARAM_SCRIND :CPX_PARAM_PREIND :CPX_PARAM_TILIM :CPX_PARAM_THREADS]
    integers= [:CPX_PARAM_SCRIND :CPX_PARAM_PREIND :CPX_PARAM_TILIM :CPX_PARAM_THREADS]
    paramDict = set_integers!(paramDict, integers)

    floats=[:CPXPARAM_MIP_Tolerances_MIPGap]
    paramDict = set_floats!(paramDict,floats)

    inputData = SolverParam(; paramDict...)

    return inputData
  catch e
    println("Can't set solver parameters.")
    throw(error())
  end
end

function set_parameters(runMode::runModeParam, case::caseData)
  try
    if runMode.setInputParameters
      println("Reading parameters from config file.")
      InputParameters = read_parameters_from_config_file("configParameters.in")
    else
      println("Reading inputparameters from JLD file.")
      InputFromFile = read_input_from_jld(case.InputPath, case.InputCase, "InputParam.jld")
      InputParameters = InputFromFile["InputParameters"]
    end

    return InputParameters
  catch e
    println("Can't set input parameters.")
    throw(error())
  end
end

function set_solverParameters()
  solverParameters = read_solverParameters_from_file()
  return solverParameters
end

# run mode (Bools)
function read_runMode_file(file = "runMode.in")

  runModeDict = read_type_to_dict(file, Bool)
  runMode = runModeParam(; runModeDict...)

  #runMode = update_default_settings!(runMode)

  println(runMode)

  return runMode
end



# hydropower system
function set_hydropower_system(runMode::runModeParam, case)
  try
    if runMode.hydroSystemFromFile
      println("Reading hydrosystem from input files")
      HY = ReadHydroData(joinpath(case.DataPath, "HydroSyst.dat"))
    else
      println("Reading hydrosystem from JLD file.")
      InputFromFile = read_input_from_jld(case.InputPath, case.InputCase, "HydroData.jld")
      HY = InputFromFile["HY"]
    end
    return HY
  catch e
    println("Can't set hydro power parameters.")
    throw(error())
  end
end

function ReadHydroData(path)
  f = open(path)
  println("\nopening file: ", path)
  for iLine = 1:3
    line = readline(f)
  end #skip first three lines
  line = readline(f)
  items = split(line, " ")
  NMod = parse(Int, items[1]) #set number of modules
  println("NMod = ", NMod)

  # Set maximum discharge per segment (DisMaxSeg) and efficiency per segment (Eff), for each hydro module
  NDSeg = zeros(Int, NMod)
  MaxSeg = 10
  DisMaxSeg = zeros(Float64, NMod, MaxSeg)
  DisPointTurb=zeros(Float64, NMod, MaxSeg)
  PowMaxSegTurb=zeros(Float64, NMod, MaxSeg)
  Eff = zeros(Float64, NMod, MaxSeg)
  m_1 = zeros(Float64, NMod)
  m_2 = zeros(Float64, NMod)
  m_3 = zeros(Float64, NMod)
  m_4 = zeros(Float64, NMod) 
  K_1_max = zeros(Float64, NMod)
  K_2_max = zeros(Float64, NMod)
  K_3_max = zeros(Float64, NMod)
  K_4_max = zeros(Float64, NMod)
  for iMod = 1:NMod
    line = readline(f)
    items = split(line, " ")
    NDSeg[iMod] = parse(Int, items[1])
    for iSeg = 1:NDSeg[iMod]
      if (iSeg == 1)
        DisMaxSeg[iMod, iSeg] = parse(Float64, items[3])
        DisPointTurb[iMod, iSeg] = parse(Float64, items[3])
        PowMaxSegTurb[iMod, iSeg] = parse(Float64, items[2])
        Eff[iMod, iSeg] = parse(Float64, items[2]) / parse(Float64, items[3])
      else
        DisMaxSeg[iMod, iSeg] = parse(Float64, items[1+iSeg*2]) - parse(Float64, items[1+(iSeg-1)*2])
        DisPointTurb[iMod, iSeg] = parse(Float64, items[1+iSeg*2])
        PowMaxSegTurb[iMod, iSeg] = parse(Float64, items[iSeg*2])
        Eff[iMod, iSeg] = (parse(Float64, items[1+iSeg*2-1]) - parse(Float64, items[1+(iSeg-1)*2-1])) / DisMaxSeg[iMod, iSeg]
      end
    end
    m_1[iMod] = (PowMaxSegTurb[iMod, 2]-PowMaxSegTurb[iMod, 1])/(DisPointTurb[iMod, 2]-DisPointTurb[iMod, 1])
    m_2[iMod] = (PowMaxSegTurb[iMod, 3]-PowMaxSegTurb[iMod, 2])/(DisPointTurb[iMod, 3]-DisPointTurb[iMod, 2])
    m_3[iMod] = (PowMaxSegTurb[iMod, 4]-PowMaxSegTurb[iMod, 3])/(DisPointTurb[iMod, 4]-DisPointTurb[iMod, 3])
    m_4[iMod] = (PowMaxSegTurb[iMod, 5]-PowMaxSegTurb[iMod, 4])/(DisPointTurb[iMod, 5]-DisPointTurb[iMod, 4])
    K_1_max[iMod] = PowMaxSegTurb[iMod, 1]-DisPointTurb[iMod, 1]*((PowMaxSegTurb[iMod, 2]-PowMaxSegTurb[iMod, 1])/(DisPointTurb[iMod, 2]-DisPointTurb[iMod, 1]))
    K_2_max[iMod] = PowMaxSegTurb[iMod, 2]-DisPointTurb[iMod, 2]*((PowMaxSegTurb[iMod, 3]-PowMaxSegTurb[iMod, 2])/(DisPointTurb[iMod, 3]-DisPointTurb[iMod, 2]))
    K_3_max[iMod] = PowMaxSegTurb[iMod, 3]-DisPointTurb[iMod, 3]*((PowMaxSegTurb[iMod, 4]-PowMaxSegTurb[iMod, 3])/(DisPointTurb[iMod, 4]-DisPointTurb[iMod, 3]))
    K_4_max[iMod] = PowMaxSegTurb[iMod, 4]-DisPointTurb[iMod, 4]*((PowMaxSegTurb[iMod, 5]-PowMaxSegTurb[iMod, 4])/(DisPointTurb[iMod, 5]-DisPointTurb[iMod, 4]))
  end

  # Set max reservoir (MaxRes), start reservoir (Resinit0), Scale and NUp for each module
  for iLine = 1:4
    line = readline(f)
  end #skip 4 lines
  MaxRes = zeros(Float64, NMod)
  ResInit0 = zeros(Float64, NMod)
  Scale = zeros(Float64, NMod)
  NUp = zeros(Int, NMod)
  qMin = zeros(Float64, NMod)
  for iMod = 1:NMod
    line = readline(f)
    items = split(line, " ")
    MaxRes[iMod] = parse(Float64, items[1])
    ResInit0[iMod] = parse(Float64, items[2])
    Scale[iMod] = parse(Float64, items[4])
    NUp[iMod] = parse(Int, items[5])
    qMin[iMod] = parse(Float64, items[6])
  end

  for iLine = 1:4
    line = readline(f)
  end #skip 4 lines

  #On which reservoir is the pump
  line=readline(f)
  items=split(line," ")
  Station_with_pump=parse(Int,items[1])

  #Creation of a pump
  NDSegPump=zeros(Int)
  MaxSegPump=10
  DisMaxSegPump=zeros(Float64,MaxSegPump)
  DisPointPump=zeros(Float64,MaxSegPump)
  PowMaxSegPump=zeros(Float64,MaxSegPump)
  EffPump=zeros(Float64,MaxSegPump)
  m_pump = 0
  K_pump_max = 0
  Pump_direction=zeros(Int,2)

  line=readline(f)
  items=split(line, " ")

  if Station_with_pump==0
    NDSegPump=0
    MaxSegPump=1
    DisMaxSegPump=zeros(Float64,MaxSegPump)
    DisPointPump=zeros(Float64,MaxSegPump)
    PowMaxSegPump=zeros(Float64,MaxSegPump)
    EffPump=zeros(Float64,MaxSegPump)
    Pump_direction=zeros(Int,2)
  elseif Station_with_pump==1
    NDSegPump=parse(Int,items[1])
    for iSeg = 1:NDSegPump
      if (iSeg == 1)
        DisMaxSegPump[iSeg] = parse(Float64, items[2])
        DisPointPump[iSeg] = parse(Float64, items[2])
        PowMaxSegPump[iSeg] = parse(Float64, items[3])
        EffPump[iSeg] = parse(Float64, items[2]) / parse(Float64, items[3])
      else
        DisMaxSegPump[iSeg] = parse(Float64, items[iSeg*2]) - parse(Float64, items[(iSeg-1)*2])
        DisPointPump[iSeg] = parse(Float64, items[iSeg*2])
        PowMaxSegPump[iSeg] = parse(Float64, items[iSeg*2+1])
        EffPump[iSeg] = DisMaxSegPump[iSeg]/(parse(Float64, items[iSeg*2+1]) - parse(Float64, items[(iSeg-1)*2+1])) 
      end
    end
    m_pump = (PowMaxSegPump[2]-PowMaxSegPump[1])/(DisPointPump[2]-DisPointPump[1])
    K_pump_max = PowMaxSegPump[1]-DisPointPump[1]*((PowMaxSegPump[2]-PowMaxSegPump[1])/(DisPointPump[2]-DisPointPump[1]))
    Pump_direction[1]=1
    Pump_direction[2]=-1
  end
  
  if Station_with_pump==0
    println("Traditional system")
  else
    println("Pumped storage system")
  end

  for iLine = 1:2
    line = readline(f)
  end #skip 3 lines
  N_min_flows=zeros(Int,NMod)
  Activation_weeks=zeros(Int,NMod,10)
  Min_flows=zeros(Float64,NMod,10)

 if any(x->x!=0,qMin)  
   for iMod=1:NMod
     line = readline(f)
     items = split(line, " ")
     N_min_flows[iMod] = parse(Int, items[1])  
      if N_min_flows[iMod] !=0
        for i=1:N_min_flows[iMod]
           Activation_weeks[iMod,i]=parse(Int,items[i+1])
        end
      end
   end

   for iLine = 1:2
     line = readline(f)
   end #skip 3 lines
  
   for iMod=1:NMod
     line = readline(f)
     items = split(line, " ")
     if N_min_flows[iMod] !=0
       for i=1:N_min_flows[iMod]
           Min_flows[iMod,i]=parse(Float64,items[i])
       end
     end
   end
 end

  close(f)
  println("closed file: ", path)

  return HydroData(NMod, NUp, Eff, NDSeg, DisMaxSeg, MaxRes, Scale, ResInit0, qMin,Station_with_pump,NDSegPump,DisMaxSegPump,DisPointPump,PowMaxSegPump,DisPointTurb,PowMaxSegTurb,EffPump,Pump_direction,N_min_flows,Activation_weeks,Min_flows,m_1,m_2,m_3,m_4,K_1_max,K_2_max,K_3_max,K_4_max,m_pump,K_pump_max)
end

# Environmental constraint
function readMaxDischargeInput(path)
  f = open(path)
  println("\nopening file: ", path)
  readline(f) #skip first line
  line = readline(f)
  items = split(line, " ")
  NEnvCon = parse(Int, items[1]) #number of constraints

  EnvCon = Vector{maxDishargeConstraint}(undef, NEnvCon)

  for con = 1:NEnvCon
    for iLine = 1:2
      line = readline(f)
    end #skip 2 lines

    line = readline(f)
    items = split(line, " ")
    envMod = parse(Int, items[1]) #set module
    firstAct = parse(Int, items[2])
    lastAct = parse(Int, items[3])
    lastMaxDisch = parse(Int, items[4])
    lastNoDecrease = parse(Int, items[5])

    for iLine = 1:2
      line = readline(f)
    end #skip 2 lines
    line = readline(f)
    items = split(line, " ")
    actLevel = parse(Float64, items[1])
    deactLevel = parse(Float64, items[2])
    maxDischarge = parse(Float64, items[3])

    EnvCon[con] = maxDishargeConstraint(
      envMod,
      firstAct,
      lastAct,
      lastMaxDisch,
      lastNoDecrease,
      actLevel,
      deactLevel,
      maxDischarge,
    )

  end

  for iLine = 1:3
    line = readline(f)
  end #skip 3 lines

  line = readline(f)
  items = split(line, " ")
  flowDependentqMin = Bool(parse(Float64, items[1]))
  if flowDependentqMin
    qMod = parse(Float64, items[2])
    qMin = parse(Float64, items[3])
    for iLine = 1:2
      line = readline(f)
    end #skip 2 lines
    line = readline(f)

    flowMin = parse.(Float64, split(line, " "))

    for iLine = 1:2
      line = readline(f)
    end #skip 2 lines
    line = readline(f)
    flowScale = parse.(Float64, split(line, " "))

    qMinDependent = qMinDependentnFlow(qMin, qMod, flowMin, flowScale)
  else
    qMinDependent = 0
  end
  return EnvCon, qMinDependent
end

function use_early_activation(envDataList)
  if first(envDataList).firstAct > 0
    return true
  else
    return false
  end
end

# case settings
function set_runCase(file = "runCase.in")

  runCaseDict = read_type_to_dict(file, String)
  case = caseData(; runCaseDict...)

  println(caseData)

  return case
end

function readPrice(PriceFile)
  f = open(PriceFile, "r")

  line = readline(f)
  line = readline(f)
  items = split(line, " ")
  NPeriods = parse(Int, items[1])
  line = readline(f)
  PriceScale = zeros(Float64, NPeriods)
  line = readline(f)
  items = split(line, " ")
  for iPeriod = 1:NPeriods
    PriceScale[iPeriod] = parse(Float64, items[iPeriod])
  end

  close(f)
  #return Price, PriceScale
  return PriceScale
end

function set_priceVar(runMode::runModeParam, case)
  #Price data
  try
    PriceScale =
      readPrice(joinpath(case.DataPath, string("prices_", case.PriceVar, ".dat")))
    return PriceScale
  catch e
    println("Can't set price.")
    throw(error())
  end
end

function read_input_from_jld(InputPath, InputCase, filname)
  inputFile = InputCase * filname
  try
    InputFromFile = load(joinpath(InputPath, inputFile))
    return InputFromFile
  catch e
    println(string("Could not read input: ", inputFile))
  end
end

# define Markov model
function set_stochastic_variables(
  runMode::runModeParam,
  case::caseData,
  InputParameters::InputParam,
)
  #Set stochastic variables, Markow model 
  try
    if runMode.createMarkovModel
      println("Creating Markov model.")
      inflow = read_csv("inflow_ita.csv", case.DataPath) #Mm3     # file degli inflow futuri
      price = read_csv("price_week.csv", case.DataPath) #Eur/MWh     # file dei prezzi futuri
      scenLData = samplingAlg(inflow, price, InputParameters)
    else
      println("Reading makov model from JLD file.")
      scenFile = "scenarioLatticeData.jld"
      ScenStatesFromFile = load(joinpath(case.InputPath, scenFile))
      scenLData = ScenStatesFromFile["scenarioLatticeData"]
    end

    return scenLData
  catch e
    println("Can't set stochastic variables, Markov model.")
    throw(error())
  end
end

function adjust_scenLattice_to_envConstriant(envDataList, runMode, scenLData)

  if !(qMinDependent == 0)
    println("State-dependent min flow included..")
  end
  println("Solving with environmental constraint..")

  if use_early_activation(envDataList)
    println("Early activation of environmental constraint is possible..")
    if runMode.createMarkovModel | runMode.extendScenarioLattice
      extendedLattice = extend_scenario_lattice(scenLData, envDataList[1])
      scenLData.ScenarioLattice.states = extendedLattice.states
      scenLData.ScenarioLattice.probability = extendedLattice.probability
    end
  end

  return scenLData
end


# define scenarios for end-simualtion
function set_sim_scenarios(runMode, NSimScen, scenLData)
  try
    if runMode.drawScenarios
      SimScen =
        drawScenForSim(NSamples, NSimScen, scenLData.trajectories, scenLData.KmeansClusters)
      println("Draw scenarios for simulation")
    end

    if runMode.drawOutofSampleScen
      SimScen = drawOutOfSampleScenForSim(priceYear, NSimScen, scenLData.states, 3456) #seednum=3456
      println("Draw out-of-sample scenarios for simulation")
    end

    if runMode.useScenariosFromDataStorage
      scenFile = "SimScen.jld"
      SimScenFromFile = load(joinpath(case.InputPath, scenFile))
      SimScen = SimScenFromFile["SimScen"]
      println("Use stored scenarios for simulation")
    end

    if runMode.useHistoricScen
      SimScen, NSimScen = scenFromInput(inflow, price, scenLData.states)
      println("Use historical scenarios for simulation")
    end

    return SimScen
  catch e
    println("Could not set simulations scenarios.")
    throw(error())
  end
end

# read results from SDP 
# -> for stand-alone sim or warmstart of SDP
function read_SDP_results(InputCase)
  try
    SDPfile = InputCase * "sdpRes.jld"
    resultsSDPfile = load(joinpath(case.InputPath, SDPfile))
    ResultsSDP = resultsSDPfile["sdpRes"]
    return ResultsSDP
  catch e
    println("Can't read future value table.")
    throw(error())
  end
end

function needToExpandAlphaTable(ResultsSDP, NStates, envDataList)
  return first(envDataList).firstAct > 0 && size(ResultsSDP.AlphaTable)[2] == NStates
end

# set case run name for saving results
function set_run_name(runMode, case, ResultPath, InputParameters)

  if runMode.envConst
    envCondition = "_wEnv"
  else
    envCondition = "_woEnv"
  end

  if runMode.ramping_constraints
    ramp = "_wRamp"
  else
    ramp = "_woRamp"
  end

  param = string("seg", InputParameters.NSeg, "_NStates", InputParameters.NStates)

  Finalpath = joinpath(ResultPath, splitdir(case.DataPath)[end])
  Finalpath = joinpath(Finalpath, param)

  if !runMode.solveMIP
    ext = "_LP"
  else
    ext = ""
  end

  runName = string(date, envCondition, ramp, case.CaseName)

  println(string("Run name: ", runName))
  println(string("Result directory: ", joinpath(Finalpath, runName)))

  return joinpath(Finalpath, runName)
end



