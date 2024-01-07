# ERROR EVALUATION

using DataFrames,XLSX

function error_evaluation(InputParameters::InputParam)

    NStep = InputParameters.NStep
    NSeg = InputParameters.NSeg
    

    Error=zeros(HY.NMod,NStage);             # for upper segment
    Upper_bound=zeros(HY.NMod,NStage);
    Lower_bound=zeros(HY.NMod,NStage);
    seg_upper= zeros(NSeg[1]);            # zeros(HY.NMod,NSeg[iMo])
    seg_lower= zeros(NSeg[2]);

   
    for iSeg=1:NSeg[1]
        seg_upper[iSeg]=(iSeg-1)* HY.MaxRes[1]/(NSeg[1]-1)
    end
    for iSeg=1:NSeg[2]
        seg_lower[iSeg]=(iSeg-1)* HY.MaxRes[2]/(NSeg[2]-1)
    end
    

    # UPPER RESERVOIR
    for iStage=1:NStage
        for iSeg=1:NSeg[1]-1                               #NSeg[iMod]
            if ResultsSim.Reservoir[1,1,iStage,NStep] >= seg_upper[iSeg] && ResultsSim.Reservoir[1,1,iStage,NStep] <= seg_upper[iSeg+1]
                    
                Upper_bound[1,iStage]= seg_upper[iSeg+1]
                Lower_bound[1,iStage]= seg_upper[iSeg]

                difference1 = Upper_bound[1,iStage]-ResultsSim.Reservoir[1,1,iStage,NStep] 
                difference2 = ResultsSim.Reservoir[1,1,iStage,NStep] - Lower_bound[1,iStage]
                    #difference2 = abs(Lower_bound[iMod,iStage]-ResultsSim.Reservoir[iMod,1,iStage,NStep])

                if difference1<= difference2
                    Error[1,iStage] = (Upper_bound[1,iStage]-ResultsSim.Reservoir[1,1,iStage,NStep])/Upper_bound[1,iStage]*100
                elseif difference1 > difference2
                    Error[1,iStage] = (ResultsSim.Reservoir[1,1,iStage,NStep]-Lower_bound[1,iStage])/ResultsSim.Reservoir[1,1,iStage,NStep]*100
                        #Error[iMod,iStage] =  abs(Lower_bound[iMod,iStage]-ResultsSim.Reservoir[iMod,1,iStage,NStep])/Lower_bound[iMod,iStage]*100
                end

            end
        end
    end

    # LOWER RESERVOIR
    for iStage=1:NStage
        for iSeg=1:NSeg[2]-1                               #NSeg[iMod]
            if ResultsSim.Reservoir[2,1,iStage,NStep] >= seg_lower[iSeg] && ResultsSim.Reservoir[2,1,iStage,NStep] <= seg_lower[iSeg+1]
                    
                Upper_bound[2,iStage]=seg_lower[iSeg+1]
                Lower_bound[2,iStage]= seg_lower[iSeg]

                difference1 = Upper_bound[2,iStage]-ResultsSim.Reservoir[2,1,iStage,NStep] 
                difference2 = ResultsSim.Reservoir[2,1,iStage,NStep] - Lower_bound[2,iStage]
                    #difference2 = abs(Lower_bound[iMod,iStage]-ResultsSim.Reservoir[iMod,1,iStage,NStep])

                if difference1<= difference2
                    Error[2,iStage] = (Upper_bound[2,iStage]-ResultsSim.Reservoir[2,1,iStage,NStep])/Upper_bound[2,iStage]*100
                elseif difference1 > difference2
                    Error[2,iStage] = (ResultsSim.Reservoir[2,1,iStage,NStep]-Lower_bound[2,iStage])/ResultsSim.Reservoir[2,1,iStage,NStep]*100
                        #Error[iMod,iStage] =  abs(Lower_bound[iMod,iStage]-ResultsSim.Reservoir[iMod,1,iStage,NStep])/Lower_bound[iMod,iStage]*100
                end

            end
        end
    end

    Upper_reservoir = DataFrame();
    Lower_reservoir = DataFrame();    

    Upper_reservoir[!,"Stages"]=1:1:NStage
    Upper_reservoir[!,"Volume last step"] = ResultsSim.Reservoir[1,1,:,end];
    Upper_reservoir[!,"Error from nearest"] = Error[1,:];
    Upper_reservoir[!,"Upper bound"] = Upper_bound[1,:];
    Upper_reservoir[!,"Lower bound"] = Lower_bound[1,:];

    Lower_reservoir[!,"Stages"]=1:1:NStage;
    Lower_reservoir[!,"Volume Last step"] = ResultsSim.Reservoir[2,1,:,end];
    Lower_reservoir[!,"Error from nearest"] = Error[2,:];
    Lower_reservoir[!,"Upper bound"] = Upper_bound[2,:];
    Lower_reservoir[!,"Lower bound"] = Lower_bound[2,:];

    

    XLSX.writetable("Error$NSeg.xlsx",overwrite=true,
    Upper_volume=(collect(DataFrames.eachcol(Upper_reservoir)),DataFrames.names(Upper_reservoir)),
    Lower_volume=(collect(DataFrames.eachcol(Lower_reservoir)),DataFrames.names(Lower_reservoir))
    )

end

