#=using ScenTrees
using Distributions
using PyPlot
=#

struct StochasticData
    State::Any
    Prob::Any
end

function addLastStagesFirst(InputData)
    #add last 2 columns first
    InputData = hcat(InputData[:, end-1:end], InputData)
    return InputData
end

function makeLattice(InputData, NlatticeStates)
    LatticeShape = repeat([NlatticeStates], NStage + 1)
    pushfirst!(LatticeShape, 1)
    Scenarios = kernel_scenarios(InputData, Normal; Markovian = true)
    Lattice = lattice_approximation(LatticeShape, Scenarios, 10000, 2)
    return Lattice
end

function getData(lattice)

    State = zeros(NStage, InfScen)
    Prob = zeros(NStage, InfScen, InfScen)
    for stage = 3:NStage+2
        for toScen = 1:size(lattice.probability[stage])[2]
            State[stage-2, toScen] = lattice.state[stage][toScen]
            for fromScen = 1:size(lattice.probability[stage])[1]
                Prob[stage-2, fromScen, toScen] =
                    lattice.probability[stage][fromScen, toScen, 1] /
                    sum(lattice.probability[stage][fromScen, :, 1])
            end
        end
    end

    return StochasticData(State, Prob)
end

function scaleInflow(lattice)

    State = zeros(HY.NMod, NStage, InfScen)
    for iMod = 1:HY.NMod
        for iStage = 1:NStage
            State[iMod, iStage, :] = lattice.State[iStage, :] .* HY.Scale[iMod]
        end
    end

    return StochasticData(State, lattice.Prob)
end

function inflowLattice(Inflow, NlatticeStates = 3, plotting = true)

    Inflow = addLastStagesFirst(Inflow)
    Lattice = makeLattice(Inflow, NlatticeStates)

    if plotting
        pygui(true)
        plt = plot_lattice(Lattice)
    end

    inflowLattice = getData(Lattice)
    inflowLattice = scaleInflow(inflowLattice)

    return inflowLattice
end



function priceLattice(Price, NlatticeStates = 3, plotting = true)

    Price = addLastStagesFirst(Price)
    Lattice = makeLattice(Price, NlatticeStates)

    if plotting
        pygui(true)
        plt = plot_lattice(Lattice)
    end

    priceLattice = getData(Lattice)

    return priceLattice
end
