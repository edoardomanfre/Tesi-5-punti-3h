# READ THE DATA FROM TXT FILE #104.10000000000001

function water_levels_evaluation(case::caseData, InputParameters::InputParam) #ResultsSim

    path=case.DataPath
    cd(path)
    f=open("Water_volumes_levels.dat")
    line=readline(f)

    line = readline(f)
    items = split(line, " ")
    NMod = parse(Int, items[1]) #set number of modules
    println("NMod = ", NMod)
    water_volumes_file=zeros(Float64,HY.NMod,21);
    water_levels_file=zeros(Float64,HY.NMod,21);
    NVolumes=zeros(NMod);
    
    for iMod=1:NMod
        line = readline(f)
        items = split(line, " ")
        NVolumes[iMod] = parse(Float64, items[1])
        for n=1:Int(NVolumes[iMod])
        water_volumes_file[iMod,n]=parse(Float64,items[1+n])
        end
    end
    water_volumes_file;

    for iLine = 1:2
        line = readline(f)
    end

    for iMod=1:NMod
        line = readline(f)
        items = split(line, " ")
        for n=1:Int(NVolumes[iMod])
        water_levels_file[iMod,n]=parse(Float64,items[n])
        end
    end
    water_levels_file;

    #EVALUATE THE WATER LEVELS, GIVEN THE WATER VOLUMES IN THE RESERVOIR
    Water_levels_Simulation=zeros(HY.NMod,NSimScen,NStage,InputParameters.NStep);
    Initial_water_levels=zeros(HY.NMod)

    # CALCULATES THE WATER LEVELS (m a.s.l) FROM THE VOLUME RESULTS
    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            for iStage=1:NStage
                for iStep=1:InputParameters.NStep

                    for n=1:Int(NVolumes[iMod])-1
                        
                        if ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] == water_volumes_file[iMod,n]
                            Water_levels_Simulation[iMod,iScen,iStage,iStep] = water_levels_file[iMod,n]
                        elseif ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep]> water_volumes_file[iMod,n] && ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep]< water_volumes_file[iMod,n+1]
                            Water_levels_Simulation[iMod,iScen,iStage,iStep] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                        end

                        if HY.ResInit0[iMod] > water_volumes_file[iMod,n] && HY.ResInit0[iMod]< water_volumes_file[iMod,n+1]
                            Initial_water_levels[iMod] =(water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(HY.ResInit0[iMod]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                        end

                    end
                    
                    if ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                        Water_levels_Simulation[iMod,iScen,iStage,iStep] = water_levels_file[iMod,Int(NVolumes[iMod])]
                    end

                end
            end
        end
    end

    #  WATER LEVEL VARIATIONS PER TIME STEP
    water_level_variations=zeros(HY.NMod,NSimScen,NStage,InputParameters.NStep);
    volume_variations=zeros(HY.NMod,NSimScen,NStage,InputParameters.NStep);
    max_min_median=zeros(HY.NMod,NSimScen,NStage,3);
    weekly_water_variations=zeros(HY.NMod,NSimScen,NStage);
    
    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            for iStage=1:NStage
                for iStep=1:InputParameters.NStep

                    if iScen==1  
                        if iStage==1 #Per la settimana t=1
                            if iStep==1  
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Initial_water_levels[iMod]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - HY.ResInit0[iMod]
                            else    #iStep>1
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Water_levels_Simulation[iMod,iScen,iStage,iStep-1]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep-1]
                            end
                        else # Per le settimane t>1
                            if iStep==1
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Water_levels_Simulation[iMod,iScen,iStage-1,end]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen,iStage-1,end]
                            else
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Water_levels_Simulation[iMod,iScen,iStage,iStep-1]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep-1]
                            end
                        end
                    else #iScen>1
                        if iStage==1
                            if iStep==1
                                water_level_variations[iMod,iScen,iStage,iStep]=Water_levels_Simulation[iMod,iScen,iStage,iStep] - Water_levels_Simulation[iMod,iScen-1,end,end]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen-1,end,end]
                            else
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Water_levels_Simulation[iMod,iScen,iStage,iStep-1]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep-1]
                            end
                        else #iStage > 1
                            if iStep==1
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Water_levels_Simulation[iMod,iScen,iStage-1,end]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen,iStage-1,end]
                            else #iStep > 1
                                water_level_variations[iMod,iScen,iStage,iStep]= Water_levels_Simulation[iMod,iScen,iStage,iStep]-Water_levels_Simulation[iMod,iScen,iStage,iStep-1]
                                volume_variations[iMod,iScen,iStage,iStep] = ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep] - ResultsSim.Reservoir_round[iMod,iScen,iStage,iStep-1]
                            end
                        end
                    end
                

                    if water_level_variations[iMod,iScen,iStage,iStep] > max_min_median[iMod,iScen,iStage,1]
                        max_min_median[iMod,iScen,iStage,1] = water_level_variations[iMod,iScen,iStage,iStep]
                    elseif water_level_variations[iMod,iScen,iStage,iStep] < max_min_median[iMod,iScen,iStage,2]
                        max_min_median[iMod,iScen,iStage,2] = water_level_variations[iMod,iScen,iStage,iStep]
                    end

                end

                max_min_median[iMod,iScen,iStage,3]=median(water_level_variations[iMod,iScen,iStage,:])
                weekly_water_variations[iMod,iScen,iStage] =sum(water_level_variations[iMod,iScen,iStage,:])
            end
        end
    end
    
    frequency=zeros(HY.NMod,NSimScen,11);

    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            frequency[iMod,iScen,1]=count(i->(i<-0.20),water_level_variations[iMod,iScen,:,:])                  # -0.2
            frequency[iMod,iScen,2]=count(i->(i>=-0.20 && i<-0.15),water_level_variations[iMod,iScen,:,:])      # -0.2 -0.15
            frequency[iMod,iScen,3]=count(i->(i>=-0.15 && i<-0.10),water_level_variations[iMod,iScen,:,:]) 
            frequency[iMod,iScen,4]=count(i->(i>=-0.10 && i<-0.05),water_level_variations[iMod,iScen,:,:])
            frequency[iMod,iScen,5]=count(i->(i>=-0.05 && i<0.00),water_level_variations[iMod,iScen,:,:]) 
            frequency[iMod,iScen,6]=count(i->(i==0),water_level_variations[iMod,iScen,:,:])  
            frequency[iMod,iScen,7]=count(i->(i>0 && i<=0.05),water_level_variations[iMod,iScen,:,:]) 
            frequency[iMod,iScen,8]=count(i->(i>0.05 && i<=0.10),water_level_variations[iMod,iScen,:,:]) 
            frequency[iMod,iScen,9]=count(i->(i>0.10 && i<=0.15),water_level_variations[iMod,iScen,:,:]) 
            frequency[iMod,iScen,10]=count(i->(i>0.15 && i<=0.20),water_level_variations[iMod,iScen,:,:]) 
            frequency[iMod,iScen,11]=count(i->(i>0.20),water_level_variations[iMod,iScen,:,:])               # 0.2
        end
    end

    return Water_levels(water_volumes_file,water_levels_file,NVolumes,Water_levels_Simulation,water_level_variations,volume_variations,max_min_median,weekly_water_variations,frequency)

end

#water_volumes_file,water_levels_file,NVolumes,