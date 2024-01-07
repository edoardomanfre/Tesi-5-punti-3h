using JuMP, PiecewiseLinearOpt
using Interpolations

function func(x, y, m = [1 2 3; 4 5 6; 7 8 9])
    itp = interpolate(m, BSpline(Linear()))
    v = itp(x, y)
    return v
end

function func_1D(x, m = [1 2 5 7 13 40])
    itp = interpolate(transpose(m), (BSpline(Linear()), NoInterp()))
    v = itp(x, 1)
    return v
end

m = Model()
@variable(m, x)
@variable(m, y)



z = piecewiselinear(m, x, y, 1:1:3, 1:1:3, (u, v) -> func(u, v))
@objective(m, Min, z)

#konklusjon: virker som picewiselinear kan brukes inn i JuMP med funksjon som lineærinterpolerer mellom punktene som input.
#Fungerer ikke kun med matrise da det ikke gir en kontinuerlig beskrivelse innenfor det definerte mulighetsrommet.
#Kan være aktuelt å bruke picewiselinear i stede for egen SOS2 implementasjon, da dette muliggjør bruk av raskere implementasjoner for triangulering.
