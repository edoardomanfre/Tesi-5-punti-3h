

using ScenTrees
using Distributions
using PyPlot
using Parameters


function testing_scen()
    Inflow = rand(5, 8)
    Price = Inflow .* 120

    LatticeShape = repeat([3], 7)
    pushfirst!(LatticeShape, 1)
    scen1 = kernel_scenarios(Inflow, Normal; Markovian = true)
    scen2 = kernel_scenarios(Price, Normal; Markovian = true)
    return hcat(scen1(), scen2())
end

Inflow = rand(5, 8)
Price = Inflow .* 120

LatticeShape = repeat([3], 7)
pushfirst!(LatticeShape, 1)
scen1 = kernel_scenarios(Inflow, Normal; Markovian = true)
a = scen1()

Lattice = lattice_approximation(LatticeShape, scen1, 10000, 2, 1)

plotting = true
if plotting
    pygui(true)
    plt = plot_lattice(Lattice)
end
