module MyTypes

export HydroData

struct HydroData
    NMod::Any
    NUp::Any
    Eff::Any         # MW/m3/s
    NDSeg::Any
    DisMaxSeg::Any   # m3/s
    MaxRes::Any      #Mm3/s
    Scale::Any
    ResInit0::Any    #Mm3
    qMin::Any
end
end
