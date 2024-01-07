# Adding pre-defined end-vaue cuts for t=T

struct EndValuation
    NEndCut::Any
    EndCutCoef::Any
    EndCutRHS::Any
end

function ReadEndValuation(NMod)
    f = open("C:\\GitSource\\models\\data\\testmodel_oneRes\\EndValueCuts.dat", "r")
    line = readline(f)
    items = split(line, " ")
    NEndCut = parse(Int, items[1])
    println("NEndValueCuts = ", NEndCut)

    EndCutCoef = zeros(NMod, NEndCut)
    EndCutRHS = zeros(NEndCut)
    for iCut = 1:NEndCut
        line = readline(f)
        items = split(line, " ")
        for iMod = 1:NMod
            EndCutCoef[iMod, iCut] = parse(Float64, items[iMod])
        end
        EndCutRHS[iCut] = parse(Float64, items[NMod+1])
    end
    close(f)
    return EndValuation(NEndCut, EndCutCoef, EndCutRHS)
end
