

using Interpolations

function linear_iterpolate(NSeg, HY, Alphatable, res_value)
    ResEnd = 0:HY.MaxRes[1]/(NSeg-1):HY.MaxRes[1]
    itp = interpolate(Alphatable, BSpline(Linear()))
    sitp = Interpolations.scale(itp, ResEnd)
    return sitp(res_value)
end
