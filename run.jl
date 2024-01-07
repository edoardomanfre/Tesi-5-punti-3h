using JuMP
using Printf
using CPLEX
using MathOptInterface
using JLD
using TimerOutputs
using Distributions
using DataFrames
using Parameters
using Dates
using CSV
using Clustering
using Random

#using PyPlot
#using FileIO
 
include("MyStructs.jl")
include("setInputParameters.jl")
include("stageprob.jl")
include("stageprob_sim.jl")
include("waterValueFunctions.jl")
include("getResults.jl")
#include("readInput.jl")
include("preprocessing.jl")
include("scenarioSampling.jl")
include("activateMaxDischargeConstraint.jl")
include("SDP.jl")
include("sim.jl")
include("PumpConstraints.jl")

include("ProductionFactors.jl")
include("WaterLevels.jl")
include("ExcelSavings.jl")
include("volume_plot.jl")
include("head_function.jl")
#include("Plots.jl")
#include("ErrorEvaluation.jl")

#path to input files and to location to save results

date = string(today())
ScalePriceMeanYear = false #if making markov model or changing markov model
 

# PREAPARE INPUT 
#---------------------------------------------------------------------------------------------
to = TimerOutput()
##@timeit to "Solve complete model" begin
@timeit to "Set input data" begin

  # Set run case
  case = set_runCase()
  @unpack (DataPath, InputPath, ResultPath, InputCase, PriceYear, PriceVar) = case

  # Set run mode
  runMode = read_runMode_file()

  InputParameters = set_parameters(runMode, case)
  @unpack (NSeg, NStage, NStates, NSamples, NSimScen, LimitPump) = InputParameters

  SolverParameters = set_solverParameters()

  FinalResPath = set_run_name(runMode, case, ResultPath, InputParameters)

  # set hydropower system
  HY = set_hydropower_system(runMode::runModeParam, case)

  if HY.NMod == 1
    ResSeg = [(s - 1) * HY.MaxRes[iMod] / (NSeg[iMod] - 1) for iMod in HY.NMod for s = 1:NSeg[iMod]]        #NSeg[iMod]
  elseif HY.NMod == 2
    AllSeg, ResSeg = twoResStateMatrix(NSeg, HY)
  else
    print("Segments can not be calculated for more than two reservoirs.")
  end  

  #inflow = set_inflow(runMode, InputFromFile, case)
  #price, PriceScale = set_price(runMode, InputFromFile, case)
  scenLData = set_stochastic_variables(runMode, case, InputParameters)
  PriceScale = read_csv("Price_ita_3h.csv",case.DataPath)                                       # price_variations_24h.csv -> vanno i valori di scala

  if ScalePriceMeanYear
    println("Changing price series mean..")
    scenLData = scalePrice(
      joinpath(case.DataPath, string("meanPrice_", priceYear, ".csv")),
      scenLData,
    )
  end

  # Set environmental constraint
  if runMode.envConst
    envDataList, qMinDependent =
      readMaxDischargeInput(joinpath(case.DataPath, "EnvConditions.dat"))
    scenLData = adjust_scenLattice_to_envConstriant(envDataList, runMode, scenLData)
  else
    qMinDependent = 0
    envDataList = []
  end

  if runMode.useWaterValues
    warmstart = load(joinpath(case.InputPath, "sdpRes.jld"))["sdpRes"].WVTable
  else
    warmstart = 0
  end
end # timer "Prep input"

# save input data 
@timeit to "Save input" begin
  save(joinpath(FinalResPath, "caseDetails.jld"), "case", case)
  save(joinpath(FinalResPath, "runMode.jld"), "runMode", runMode)
  save(joinpath(FinalResPath, "SolverParam.jld"), "SolverParameters", SolverParameters)
  save(joinpath(FinalResPath, "InputParam.jld"), "InputParameters", InputParameters)
  save(joinpath(FinalResPath, "HydroData.jld"), "HY", HY)
  save(joinpath(FinalResPath, "PriceScale.jld"), "PriceScale", PriceScale)
  save(joinpath(FinalResPath, "scenarioLatticeData.jld"), "scenarioLatticeData", scenLData)
end # timer save input

# STOCHASTIC OPTIMISATION
#---------------------------------------------------------------------------------------------

# Run SDP algorithm to generate expected future value table
if runMode.solveSDP
  @timeit to "SDP algorithm" ResultsSDP, stageprobResult = SDP(
    InputParameters,
    SolverParameters,
    HY,
    scenLData,
    ResSeg,
    PriceScale,
    envDataList,
    runMode,
    FinalResPath,
    warmstart,
  )
  # save SDP results
  save(joinpath(FinalResPath, "sdpRes.jld"), "sdpRes", ResultsSDP)
else
  "SDP algorithm not solved."
end

# DETERMINISTIC END SIMULATION
#---------------------------------------------------------------------------------------------

# Run simulation using expected future value curves from SDP
if runMode.simulate
  SimScen = set_sim_scenarios(runMode, NSimScen, scenLData)
  
  save(joinpath(FinalResPath, "SimScen.jld"), "SimScen", SimScen)
  if !runMode.solveSDP
    ResultsSDP = read_SDP_results(InputCase)
  end
  if runMode.envConst && needToExpandAlphaTable(ResultsSDP, NStates, envDataList)
    ResultsSDP = expand_AlphaTable(ResultsSDP)
  end

  @timeit to "Simulation" begin
    if runMode.parallellSim
      ResultsSim = paraSim(InputParameters, SolverParameters, ResultsSDP, SimScen, runMode)
    else
      ResultsSim = sim(InputParameters, SolverParameters, ResultsSDP, SimScen, runMode   )
    end
  end #timer "SDP simulation"  

  # save sim results
  save(joinpath(FinalResPath, "simRes.jld"), "simRes", ResultsSim)

else
  println("Solved without end-simulation.")
end

# TIME FOR WEATER LEVELS EVALUATION

if runMode.water_levels
    Results_Water_levels = water_levels_evaluation(case,InputParameters)
    save(joinpath(FinalResPath,"Water_levels.jld"), "WaterLevelsAnalysis", Results_Water_levels)
end

#TIME PRODUCTION FACTOR AND TIME EVALUATION

if runMode.production_factors
  Production_factors = production_factors(InputParameters)
  save(joinpath(FinalResPath,"Production_factors.jld"), "ProductionFactor", Production_factors)
end

# TIME TO SAVE DATA IN EXCEL FORMAT

if runMode.error_evaluation 
  
  path=case.ResultPath
  cd(path)

  error_evaluation(InputParameters)

end

if runMode.save_excel

  path=case.ResultPath
  cd(path)
  
  Save_data= data_saving(runMode,SimScen, InputParameters)

end

savePlots(InputParameters, ResultsSim) 

print(to)


  