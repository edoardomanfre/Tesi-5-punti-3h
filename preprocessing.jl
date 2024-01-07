function twoResStateMatrix(NSeg::Any, HY::HydroData)
    allSeg = collect(Iterators.product(1:NSeg[1], 1:NSeg[2]))                     # (1:NSeg[1],1:NSeg[2])   
    resState = Array{Tuple{Float64,Float64},2}(undef, NSeg[1], NSeg[2])           # (undef,NSeg[1],NSeg[2])
    #resState =  fill(Float64[], NSeg, NSeg)
    for i = 1:NSeg[1]                                                          # NSeg[1]
        for j = 1:NSeg[2]                                                      # NSeg[2]
            resState[i, j] =
                ((i - 1) * HY.MaxRes[1] / (NSeg[1] - 1), (j - 1) * HY.MaxRes[2] / (NSeg[2] - 1))          
            #resState[i,j] = [(i-1)*HY.MaxRes[1]/(NSeg[1]-1), (j-1)*HY.MaxRes[2]/(NSeg[2]-1)]
        end
    end
    return allSeg, resState
end


function expand_AlphaTable(ResultsSDP) 
    println("Expanding AlphaTable to consider early activation")
    sizeAlphaTable = size(ResultsSDP.AlphaTable)
    newAlphaTable = zeros(sizeAlphaTable[1], sizeAlphaTable[2]*2, sizeAlphaTable[3])
    for t = 1 : sizeAlphaTable[1]
        for resSeg = 1 : sizeAlphaTable[3]
            newAlphaTable[t,1:sizeAlphaTable[2],resSeg] .= ResultsSDP.AlphaTable[t,:,resSeg]
            newAlphaTable[t,sizeAlphaTable[2]+1:end, resSeg] .= ResultsSDP.AlphaTable[t,:,resSeg]
        end
    end
    return FutureValueSDP(ResultsSDP.ResSeg,  ResultsSDP.WVTable, newAlphaTable)
end