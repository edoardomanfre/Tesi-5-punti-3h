

function ReadHydroData(path)
  f = open(path)
  println("opening file: ", path)
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
  Eff = zeros(Float64, NMod, MaxSeg)
  for iMod = 1:NMod
    line = readline(f)
    items = split(line, " ")
    NDSeg[iMod] = parse(Int, items[1])
    for iSeg = 1:NDSeg[iMod]
      if (iSeg == 1)
        DisMaxSeg[iMod, iSeg] = parse(Float64, items[1+iSeg*2])
        Eff[iMod, iSeg] = parse(Float64, items[2]) / parse(Float64, items[3])
      else
        DisMaxSeg[iMod, iSeg] =
          parse(Float64, items[1+iSeg*2]) - parse(Float64, items[1+(iSeg-1)*2])
        Eff[iMod, iSeg] =
          (parse(Float64, items[1+iSeg*2-1]) - parse(Float64, items[1+(iSeg-1)*2-1])) /
          DisMaxSeg[iMod, iSeg]
      end
    end
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

  #On which reservoir the pump is installed
  line=readline(f)
  items=split(line," ")
  Station_with_pump=parse(Int,items[1])

  #Creation of a pump
  NDSegPump=zeros(Int)
  MaxSegPump=10
  DisMaxSegPump=zeros(Float64,MaxSegPump)
  EffPump=zeros(Float64,MaxSegPump)
  Pump_direction=zeros(Int,2)

  line=readline(f)
  items=split(line, " ")

  if Station_with_pump==0
    NDSegPump=0
    MaxSegPump=1
    DisMaxSegPump=zeros(Float64,MaxSegPump)
    EffPump=zeros(Float64,MaxSegPump)
    Pump_direction=zeros(Int,2)
  elseif Station_with_pump==1 
    NDSegPump=parse(Int,items[1])
    for iSeg = 1:NDSegPump
      if (iSeg == 1)
        DisMaxSegPump[iSeg] = parse(Float64, items[2])
        EffPump[iSeg] = parse(Float64, items[2]) / parse(Float64, items[3])
      else
        DisMaxSegPump[iSeg] = parse(Float64, items[iSeg*2]) - parse(Float64, items[(iSeg-1)*2])
        EffPump[iSeg] = DisMaxSegPump[iSeg]/(parse(Float64, items[iSeg*2+1]) - parse(Float64, items[(iSeg-1)*2+1])) 
      end
    end
    Pump_direction[1]=1
    Pump_direction[2]=-1
  end

  close(f)
  println("closed file: ", path)

  return HydroData(NMod, NUp, Eff, NDSeg, DisMaxSeg, MaxRes, Scale, ResInit0, qMin,Station_with_pump, NDSegPump, DisMaxSegPump,EffPump,Pump_direction,N_min_flows,min_flow)
end

function readInflow(InflowFile)
  f = open(InflowFile, "r")
  line = readline(f)
  line = readline(f)
  items = split(line, " ")
  NSeries = parse(Int, items[1])
  line = readline(f)

  Inflow = zeros(Float64, NSeries, NStage)

  for iStage = 1:NStage
    line = readline(f)
    items = split(line, " ")
    for iSeries = 1:NSeries
      Inflow[iSeries, iStage] = parse(Float64, items[iSeries+1])
    end
  end
  close(f)
  return Inflow
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

function read_input_from_results(InputPath, InputCase)
  inputFile = InputCase * "Input.jld"
  try
    InputFromFile = load(joinpath(InputPath, inputFile))
    return InputFromFile
  catch e
    println(string("Could not read input from previous results ", inputFile))
  end
end

function set_runCase(file = "runCase.in")

  runCaseDict = read_type_to_dict(file, String)
  case = caseData(; runCaseDict...)

  println(caseData)

  return case
end

function set_inflow(runMode::runModeParam, InputFromFile, case)
  # Inflow data 
  try
    if runMode.inflowFromFile && !runMode.inflowFromDataStorage
      println("Reading inflow from file.")
      inflow = Matrix{Float64}(
        CSV.read(joinpath(case.DataPath, "inflow_ita.csv"), DataFrame, header = false),
      ) #mm3     
    end

    if runMode.inflowFromDataStorage && !runMode.inflowFromFile
      println("Reading inflow from JLD file.")
      inflow = InputFromFile["inflow"]
    end

    return inflow
  catch e
    println("Can't set inflow.")
    throw(error())
  end
end

function set_price(runMode::runModeParam, InputFromFile, case)
  #Price data
  try
    if runMode.priceFromFile && !runMode.priceFromDataStorage
      println("Reading price from file.")
      price = Matrix{Float64}(
        CSV.read(joinpath(case.DataPath, "price_week.csv"), DataFrame, header = false),
      ) #Eur/MWh   
      PriceScale =
        readPrice(joinpath(case.DataPath, string("prices_", case.PriceVar, ".dat")))
    end

    if runMode.priceFromDataStorage && !runMode.priceFromFile
      println("Reading price from JLD file.")
      price = InputFromFile["price"]
      PriceScale = InputFromFile["PriceScale"]
    end

    return price, PriceScale
  catch e
    println("Can't set price.")
    throw(error())
  end
end

function set_stochastic_variables(
  runMode::runModeParam,
  InputFromFile,
  inflow,
  price,
  InputParameters,
)
  #Set stochastic variables, Markow model 
  try
    if runMode.createMarkovModel && !runMode.markovModelFromDataStorage
      println("Creating Markov model.")
      scenLData = samplingAlg(inflow, price, InputParameters)
    end

    if runMode.markovModelFromDataStorage && !runMode.createMarkovModel
      println("Reading makov model from JLD file.")
      scenLData = InputFromFile["scenarioLatticeData"]
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

function use_early_activation(envDataList)
  if first(envDataList).firstAct > 0
    return true
  else
    return false
  end
end

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
      SimScen = InputFromFile["SimScen"]
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
