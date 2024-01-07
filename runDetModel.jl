#Run deterministic model


using JuMP
using Printf
using CPLEX
using MathOptInterface
using FileIO
using JLD
using TimerOutputs
using Distributions
#using PyPlot
using DataFrames
using Parameters
using Dates
using Test

include("setInputParameters.jl")
include("stageprob.jl")
include("waterValueFunctions.jl")
#include("resultFunctions.jl")
include("readinput.jl")
include("preprocessing.jl")
include("scenarioSampling.jl")
include("deterministicModel.jl")
include("activateMaxDischargeConstraint.jl")

#Option to save results from SP
DebugSP = false

to = TimerOutput()
@timeit to "SDP 2 res model" begin
    @timeit to "Set input data" begin

        #path to input files and to location to save results

        DataPath = "C:\\GitSource\\models\\data\\Deterministisk_MIP\\enmag_test_MIP_early_simplePrice"
        InputPath = "C:\\Users\\linnesc\\OneDrive - NTNU\\Modellering\\Resultater\\testing\\DeterministiskMIP"
        Resultspath = "C:\\Users\\linnesc\\OneDrive - NTNU\\Modellering\\Resultater\\testing\\DeterministiskMIP"

        inputcase = "" #2021-04-23highpriceW26_lowInflow_Uenv.jld"#2021-03-23_1M_IK_P30_30Seg_5scen_Uenv_withoutEarlyActivation_noMinFlow.jld"#"2021-03-18_1M_IK_P30_5Seg_5scen_Uenv_histsim_withoutEarlyActivation.jld"# "2021-03-18_1M_IK_P30_5Seg_5scen_Uenv_histsim_withoutEarlyActivationv.jld"
        inputFile = "Input_" * inputcase #"Input_Bergsdalen_210226_tomag_kaldstad_2030_10seg_LP.jld"# "Input_Bergsdalen_210218_enmag.jld" #false
        resultTwoRes = "simRes_" * inputcase

        date = string(today())
        runName = date * "_early_simplePrice"
        envConst = true
        solveMIP = true
        earlyActivation = false


        if length(inputcase) < 1
            InputParameters = setInputParam()                                                                           # Seto i parametri di input
            @unpack (
                NSeg,
                NPrice,
                NStage,
                NHoursStage,
                NHoursStep,
                Big,
                AlphaMax,
                PenSpi,
                PenQ,
                NStep,
                MM3Week,
                MM3Step,
                StepFranc,
                NStates,
                MaxIt,
                conv,
                NSamples,
                NSimScen,
            ) = InputParameters

            #HY = ReadHydroData(joinpath(DataPath, "HydroSyst.dat"))

            inflow = Matrix{Float64}(
                CSV.read(joinpath(DataPath, "inflow.csv"), DataFrame, header = false),
            ) #mm3
            price = Matrix{Float64}(
                CSV.read(joinpath(DataPath, "price.csv"), DataFrame, header = false),
            ) #Eur/MWh

            NSimScen = size(inflow)[2]                                                                                      # indica numero degli scenari
            SimScen = SimScenData([hcat(inflow, price)], [])                                                                # per ogni scenario , ci sono 52 settimane e ogni settimana ha 1 solo valore di Inflow con relativo prezzo           
        else
            InputFromFile = load(joinpath(InputPath, inputFile))
            InputParameters = InputFromFile["InputParameters"]
            @unpack (
                NSeg,
                NPrice,
                NStage,
                NHoursStage,
                NHoursStep,
                Big,
                AlphaMax,
                PenSpi,
                PenQ,
                NStep,
                MM3Week,
                MM3Step,
                StepFranc,
                NStates,
                MaxIt,
                conv,
                NSamples,
                NSimScen,
            ) = InputParameters

            #HY = InputFromFile["HY"]
            inflow = InputFromFile["inflow"]
            price = InputFromFile["price"]
            PriceScale = InputFromFile["PriceScale"]
            scenLData = InputFromFile["scenarioLatticeData"]

            SimScen = InputFromFile["SimScen"]

        end

        @test NStep == 1

        HY = ReadHydroData(joinpath(DataPath, "HydroSyst.dat"))
        if envConst
            envData = readMaxDischargeInput(joinpath(DataPath, "EnvConditions.dat"))
        end

        if HY.NMod == 1
            ResSeg =
                [(s - 1) * HY.MaxRes[iMod] / (NSeg - 1) for iMod in HY.NMod for s = 1:NSeg]
        elseif HY.NMod == 2
            AllSeg, ResSeg = twoResStateMatrix(NSeg, HY)
        else
            print("It is not yet possible to have more than two reservoirs.")
        end
    end #timer "Set input data"

    #Run detmonistic model
    @timeit to "Optimising deterministic problem" begin

        ResultsDetModel = detOpt(SimScen, envConst)
    end #timer "SDP simulation"

    #Save results to Resultspath
    @timeit to "Save results" begin

        save(
            joinpath(Resultspath, string("Input_" * runName, ".jld")),
            "InputParameters",
            InputParameters,
            "HY",
            HY,
            "inflow",
            inflow,
            "price",
            price,
            "SimScen",
            SimScen,
        )

        save(
            joinpath(Resultspath, string("detRes_" * runName, ".jld")),
            "detRes",
            ResultsDetModel,
        )
    end #timer "Save results"
end#timer "SDP 2 res model"
print(to)
