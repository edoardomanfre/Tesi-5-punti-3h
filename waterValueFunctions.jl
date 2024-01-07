
function check_diff(Table)
    if size(Table)[end] == 1
        diff_matrix = Table[1, :, :, :] .- Table[end, :, :, :]
        #diff = [sum(abs.(diff_matrix))]
        diff = abs.(diff_matrix)
    elseif size(Table)[end] == 2
        diff_u = Table[1, :, :, 1] .- Table[end, :, :, 1]
        diff_l = Table[1, :, :, 2] .- Table[end, :, :, 2]
        #diff = [sum(abs.(diff_u)), sum(abs.(diff_l))]
        diff = [abs.(diff_u); abs.(diff_l)]
    end
    return diff
end
#=
function update_futureValue_oneres(AlphaTable,WVTable,ResSeg)
    # NB: Ulik utforming for vannverditabell 1 eller 2 magasiner
    #=
    #Alternativ Formulering
    AlphaTable[end,:,1] = DValTable[1, :]
    AlphaTable[end,:,2:end] = AlphaTable[1, :, 2:end]
    =#
    if WVTable[1,1,1] isa Float64
        #print("one-res model")
        AlphaTable[end,:,1] .= 0
        for iSeg =2:NSeg
            AlphaTable[end,:, iSeg] =
                (ResSeg[iSeg]-ResSeg[iSeg-1]) .*WVTable[1, :, iSeg-1] +
                AlphaTable[end,:, iSeg-1]
        end
    elseif WVTable[1,1,1] isa Tuple
        UpperValue = zeros( NStates,  NSeg, NSeg) #change in value from segment to segment
        LowerValue = zeros( NStates,  NSeg, NSeg) #change in value from segment to segment
        AlphaTable[end,:,1] .= 0
        for iSeg =2:length(ResSeg)
            segLower = Int(ceil(iSeg/NSeg))
            segUpper = iSeg - (segLower-1)*NSeg

            println("segLower:", segLower, " segUpper ", segUpper )

            #println("iSeg:", iSeg,"row:", row, "c:", column)
            for state =1:NStates
                if segLower == 1
                    UpperValue[state, segUpper, segLower] = (ResSeg[iSeg][1]-
                        ResSeg[iSeg-1][1])*WVTable[1, state, iSeg][1] +
                        UpperValue[state, segUpper-1, segLower]
                elseif segUpper == 1
                    LowerValue[state, segUpper, segLower] = (ResSeg[iSeg][2]-
                        ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2] +
                        LowerValue[state, segUpper, segLower-1]
                else
                    UpperValue[state, segUpper, segLower] = (ResSeg[iSeg][1]-
                        ResSeg[iSeg-1][1]) *WVTable[1, state, iSeg][1] +
                        UpperValue[state, segUpper-1, segLower]

                    LowerValue[state,  segUpper, segLower] = (ResSeg[iSeg][2]-
                        ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2] +
                        LowerValue[state,  segUpper, segLower-1]
                end
                AlphaTable[end, state, iSeg] =
                    UpperValue[ state,  segUpper, segLower] +
                    LowerValue[state,  segUpper, segLower]
            end
        end
        #println(AlphaTable[end, 5, :])
    end
    return  AlphaTable, UpperValue, LowerValue
end
=#
function update_futureValue(WVTable, ResSeg, NEnvStates)

    Table = zeros(NStates * NEnvStates, length(ResSeg))

    if size(WVTable)[end] == 2
        for iSeg = 2:length(ResSeg)
            segLower = Int(ceil(iSeg / NSeg)) #column number
            segUpper = iSeg - (segLower - 1) * NSeg #row number
            if segLower == 1 
                Table[:, iSeg] .=
                    Table[:, iSeg-1] .+
                    (ResSeg[iSeg][1] - ResSeg[iSeg-1][1]) .* WVTable[1, :, iSeg, 1]
                #AlphaTable_extra[end, state, iSeg] = AlphaTable_extra[end, state, iSeg-1] + (ResSeg[iSeg][1]-
                #        ResSeg[iSeg-1][1])*WVTable[1, state, iSeg][1]
            elseif segUpper == 1
                Table[:, iSeg] .=
                    Table[:, iSeg-NSeg] .+
                    (ResSeg[iSeg][2] - ResSeg[iSeg-NSeg][2]) .* WVTable[1, :, iSeg, 2]
                #AlphaTable_extra[end, state, iSeg] =  AlphaTable_extra[end, state, iSeg-NSeg] +(ResSeg[iSeg][2]-
                #        ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2]
            else
                Table[:, iSeg] .=
                    Table[:, iSeg-1] .+
                    (ResSeg[iSeg][1] - ResSeg[iSeg-1][1]) .* WVTable[1, :, iSeg, 1]
                #AlphaTable_extra[end, state, iSeg]  = AlphaTable_extra[end, state, iSeg-NSeg]  + (ResSeg[iSeg][2]-
                #    ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2]
            end
        end
    elseif size(WVTable)[end] == 1
        iMod = 1
        for iSeg = 2:length(ResSeg)
            Table[:, iSeg] .=
                Table[:, iSeg-1] .+
                (ResSeg[iSeg][iMod] - ResSeg[iSeg-1][iMod]) .* WVTable[1, :, iSeg, iMod]
        end
    end

    return Table
end

function calculateWatervalues(nfrom, ResSeg, AlphaTable)

    Table = zeros(HY.NMod)

    if HY.NMod == 2
        if nfrom == 1
            Table[:] .= 0
        elseif nfrom <= NSeg
            Table[1] =
                (AlphaTable[nfrom] - AlphaTable[nfrom-1]) /
                (ResSeg[nfrom][1] - ResSeg[nfrom-1][1])
            Table[2] = 0
        elseif ((nfrom - 1) / NSeg) % 1 == 0
            Table[1] = 0
            Table[2] =
                (AlphaTable[nfrom] - AlphaTable[nfrom-NSeg]) /
                (ResSeg[nfrom][2] - ResSeg[nfrom-NSeg][2])
        else
            Table[1] =
                (AlphaTable[nfrom] - AlphaTable[nfrom-1]) /
                (ResSeg[nfrom][1] - ResSeg[nfrom-1][1])
            Table[2] =
                (AlphaTable[nfrom] - AlphaTable[nfrom-NSeg]) /
                (ResSeg[nfrom][2] - ResSeg[nfrom-NSeg][2])
        end
    elseif HY.NMod == 1
        iMod = 1
        if nfrom == 1
            Table[iMod] = 0
        else
            Table[iMod] =
                (AlphaTable[nfrom] - AlphaTable[nfrom-1]) /
                (ResSeg[nfrom][iMod] - ResSeg[nfrom-1][iMod])
        end
    else
        error(" Number of modules not supported: ", HY.NMod)
    end
    return Table
end

function isNotConvex(WVTable)
    notConvex = false
    if HY.NMod == 1
        for iSeg = 3:length(ResSeg)
            if any(WVTable[iSeg-1, :] .< WVTable[iSeg, :])
                notConvex = true
                break
            end
        end
    elseif HY.NMod == 2
        for iSeg = 2:length(ResSeg)
            if (
                iSeg > 2 * NSeg &&
                ((iSeg - 1) / NSeg) % 1 != 0 &&
                ((iSeg - 2) / NSeg) % 1 != 0
            )
                if (
                    any(WVTable[iSeg-1, :] .< WVTable[iSeg, :]) ||
                    any(WVTable[iSeg-NSeg, :] .< WVTable[iSeg, :])
                )
                    #println("if 1.., ", iSeg)
                    notConvex = true
                    break
                end
            elseif (iSeg > NSeg && ((iSeg - 1) / NSeg) % 1 != 0)
                if (
                    any(WVTable[iSeg-1, 2] .< WVTable[iSeg, 2]) ||
                    any(WVTable[iSeg-NSeg, 1] .< WVTable[iSeg, 1])
                )
                    #println("if 1.., ", iSeg)
                    notConvex = true
                    break
                end
            elseif iSeg <= NSeg
                if (any(WVTable[iSeg-1, 2] .< WVTable[iSeg, 2]))
                    #println("if 1.., ", iSeg)
                    notConvex = true
                    break
                end

            elseif ((iSeg - 1) / NSeg) % 1 == 0
                if any(WVTable[iSeg-NSeg, 1] .< WVTable[iSeg, 1])
                    #println("if 1.., ", iSeg)
                    notConvex = true
                    break
                end
            end
        end
    end
    return notConvex
end

function isNotConvex2(WVTable)
    notConvex = false
    if HY.NMod == 1
        for iSeg = 3:length(ResSeg)
            if any(WVTable[iSeg-1, :] .< WVTable[iSeg, :])
                notConvex = true
                break
            end
        end
    elseif HY.NMod == 2
        for iSeg = 2:length(ResSeg)
            if iSeg > NSeg
                if WVTable[iSeg-NSeg, 1] != 0
                    if (any(WVTable[iSeg-NSeg, 1] .< WVTable[iSeg, 1]))
                        notConvex = true
                        break
                    end
                end

                if WVTable[iSeg-NSeg, 2] != 0
                    if (any(WVTable[iSeg-NSeg, 2] .< WVTable[iSeg, 2]))
                        notConvex = true
                        break
                    end
                end
            elseif (((iSeg - 1) / NSeg) % 1 != 0)
                if WVTable[iSeg-1, 1] != 0
                    if (any(WVTable[iSeg-1, 1] .< WVTable[iSeg, 1]))
                        notConvex = true
                        break
                    end
                end
                if WVTable[iSeg-1, 2] != 0
                    if (any(WVTable[iSeg-1, 2] .< WVTable[iSeg, 2]))
                        notConvex = true
                        break
                    end
                end
            end
        end
    end
    return notConvex
end
