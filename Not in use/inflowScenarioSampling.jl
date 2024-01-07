#using Random
#using Distributions
#using StatsBase


struct ScenTree
    Inflow::Any
    InfProb::Any
end

function InflowScen(NStage, InfScen)

    #Initialize randon number generation
    Random.seed!(123)

    #Read inflow Statistics
    InflowMean = zeros(Float64, NStage)
    InflowSDev = zeros(Float64, NStage)

    f = open("C:\\GitSource\\models\\data\\testmodel_oneRes\\InflowRecords.dat", "r")
    for iLine = 1:3
        line = readline(f)
    end
    line = readline(f)
    items = split(line, " ")
    ACorr = parse(Float64, items[1])
    line = readline(f)
    for iStage = 1:NStage
        line = readline(f)
        items = split(line, " ")
        InflowMean[iStage] = parse(Float64, items[2])
    end
    line = readline(f)
    line = readline(f)
    for iStage = 1:NStage
        line = readline(f)
        items = split(line, " ")
        InflowSDev[iStage] = parse(Float64, items[2])
    end
    close(f)

    NormInflowDelta = zeros(NStage, InfScen)
    NormInflow = zeros(NStage, InfScen)
    Inflow = zeros(NStage, InfScen)
    #infProb =  zeros(NStage, InfScen)
    for iStage = 1:NStage
        NDistInflow = Normal(0.0, InflowSDev[iStage])
        #test = rand(NDistInflow, InfScen)
        NormInflowDelta[iStage, :] = rand(NDistInflow, InfScen)
        for scen = 1:InfScen
            if iStage == 1
                NormInflow[iStage, scen] = NormInflowDelta[iStage, scen]
                Inflow[iStage, scen] =
                    NormInflow[iStage, scen] * InflowSDev[iStage] + InflowMean[iStage]
                #Inflow[iStage, scen]  = scen*4
            else
                #NormInflow[iStage, scen] = ACorr*NormInflow[iStage-1, scen] + NormInflowDelta[iStage, scen]
                NormInflow[iStage, scen] =
                    0.4 * NormInflow[iStage-1, scen] + NormInflowDelta[iStage, scen]
                Inflow[iStage, scen] =
                    NormInflow[iStage, scen] * InflowSDev[iStage] + InflowMean[iStage]
                #Inflow[iStage, scen]  = scen*4
            end
        end
    end
    Inflow[Inflow.<0] .= 0

    InfProb = zeros(InfScen, InfScen)
    for fromScen = 1:InfScen
        for toScen = 1:InfScen
            InfProb[fromScen, toScen] = 1 / InfScen
        end
    end


    return ScenTree(Inflow, InfProb)
end
