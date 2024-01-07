function index_segmentChange(nfrom, d)
    nfrom_prev = []
    for i =1:length(size(AllSeg))
        if i == d append!(nfrom_prev,nfrom[i]-1)
        else append!(nfrom_prev,nfrom[i])
        end
    end
end

function check_diff(Table)
    if Table[1,1,1] isa Float64
        Diff_matrix = Table[1, :, :] - Table[end, :, :]
        diff = [sum(abs.(diff_matrix))]
    elseif Table[1,1,1] isa Tuple
        diff_u = [Table[1, state, res][1]-Table[end,state, res][1]
            for state = 1:NStates for res = 1:length(ResSeg)]
        diff_l = [Table[1, state, res][2]-Table[end,state, res][2]
            for state = 1:NStates for res = 1:length(ResSeg)]
        diff = [sum(abs.(diff_u)), sum(abs.(diff_l))]
    end
    return diff
end

function update_futureValue(AlphaTable,WVTable,ResSeg)
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
        UpperValue = zeros( NStates,  NSeg, NSeg) #change in value from segment to segment          #NSeg*7
        LowerValue = zeros( NStates,  NSeg, NSeg) #change in value from segment to segment          #NSeg*7
        AlphaTable[end,:,1] .= 0
        for iSeg =2:length(ResSeg)
            segLower = Int(ceil(iSeg/NSeg))                         #NSeg*7
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

                    LowerValue[state,  segUpper, segLower] = (ResSeg[iSeg][2]-          #NSeg*7
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

function update_futureValue_2(AlphaTable,WVTable,ResSeg)
    if WVTable[1,1,1] isa Tuple
        #AlphaTable_extra = copy(AlphaTable)
        AlphaTable[end,:,1] .= 0
        #AlphaTable_extra[end,:,1] .= 0

        for iSeg =2:length(ResSeg)
            segLower = Int(ceil(iSeg/NSeg))                                         #NSeg*7
            segUpper = iSeg - (segLower-1)*NSeg

            #println("segLower:", segLower, " segUpper ", segUpper )

            #println("iSeg:", iSeg,"row:", row, "c:", column)
            for state =1:NStates
                if segLower == 1
                    AlphaTable[end, state, iSeg] = AlphaTable[end, state, iSeg-1] + (ResSeg[iSeg][1]-
                        ResSeg[iSeg-1][1])*WVTable[1, state, iSeg][1]
                    #AlphaTable_extra[end, state, iSeg] = AlphaTable_extra[end, state, iSeg-1] + (ResSeg[iSeg][1]-
                    #        ResSeg[iSeg-1][1])*WVTable[1, state, iSeg][1]
                elseif segUpper == 1
                    AlphaTable[end, state, iSeg] =  AlphaTable[end, state, iSeg-NSeg]+ (ResSeg[iSeg][2]-                                    #NSeg*7
                        ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2]
                    #AlphaTable_extra[end, state, iSeg] =  AlphaTable_extra[end, state, iSeg-NSeg] +(ResSeg[iSeg][2]-
                    #        ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2]
                else
                    AlphaTable[end, state, iSeg] = AlphaTable[end, state, iSeg-1] +(ResSeg[iSeg][1]-
                        ResSeg[iSeg-1][1]) *WVTable[1, state, iSeg][1]

                    #AlphaTable_extra[end, state, iSeg]  = AlphaTable_extra[end, state, iSeg-NSeg]  + (ResSeg[iSeg][2]-
                    #    ResSeg[iSeg-NSeg][2])*WVTable[1, state, iSeg][2]
                end
            end
        end
        #println(AlphaTable[end, 5, :])
    end

    #diff = AlphaTable[end,:,:] .- AlphaTable_extra[end,:,:]

    #if findmax(diff)[1] < 0.00001
    #    println("Alpha tables equal")
    #else
    #    println("Alpha tables not equal")
    #end
    #return  AlphaTable, AlphaTable_extra
    return  AlphaTable
end

function make_stageprobResult()
    stageprobResult = Dict("obj" => zeros(MaxIt, NStage, length(ResSeg), NStates),
        "inflow" => zeros(MaxIt, NStage, length(ResSeg), NStates, HY.NMod),
        "resInit" => zeros(MaxIt, NStage, length(ResSeg), NStates, HY.NMod),
        "alphanfrom" => zeros(MaxIt, NStage, length(ResSeg), NStates),
        "Reservoir"  => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "Spillage"  => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "Production"  => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "Tanking"  => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod),
        "disSeg"  => [zeros(MaxIt, NStage, length(ResSeg), NStates, NStep,
            HY.NSeg[iMod]) for iMod=1:HY.NMod],
        "gamma" => zeros(MaxIt, NStage, length(ResSeg), NStates, NSeg, NSeg),
        "totDischarge" => zeros(MaxIt, NStage, length(ResSeg), NStates, NStep, HY.NMod))
    return stageprobResult
end



function check_results(n,t, nfrom, state, stageprobResult, NMod, SP::StageProblem)
    stageprobResult["obj"][n, t, nfrom, state] = JuMP.objective_value(SP.model)
    stageprobResult["alphanfrom"][n, t, nfrom, state] = JuMP.value(SP.alpha)
    for iMod =1:NMod
        stageprobResult["inflow"][n, t, nfrom, state, iMod] =
            JuMP.value(SP.inf[iMod])
        stageprobResult["resInit"][n, t, nfrom, state, iMod] =
            JuMP.value(SP.resbalInitRHS[iMod])
        for iStep = 1:NStep
            stageprobResult["Reservoir"][n, t, nfrom, state, iStep, iMod] =
                JuMP.value(SP.res[iMod,iStep])
            stageprobResult["Spillage"][n, t,  nfrom, state, iStep, iMod] =
                JuMP.value(SP.spill[iMod,iStep])
            stageprobResult["Production"][n, t, nfrom, state, iStep,iMod] =
                JuMP.value(SP.prod[iMod,iStep])
            stageprobResult["Tanking"][n, t, nfrom, state, iStep, iMod] =
                JuMP.value(SP.tank[iMod,iStep])
            for iSeg = 1:HY.NSeg[1]
                stageprobResult["disSeg"][iMod][n, t, nfrom, state, iStep, iSeg] =
                    JuMP.value(SP.disSeg[iMod, iSeg, iStep])
            end
            stageprobResult["totDischarge"][n, t, nfrom, state, iStep,iMod] =
                sum(stageprobResult["disSeg"][iMod][n, t, nfrom, state, iStep, :])
        end
    end

    for nU =1:NSeg
        for nL=1:NSeg                   #NSeg*7
            stageprobResult["gamma"][n, t, nfrom, state, nU,nL] =
                JuMP.value(SP.gamma[nU,nL])
        end
    end
    return stageprobResult
end
