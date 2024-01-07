# Function used to calculate the head

function head_evaluation(
    case::caseData, 
    Reservoir_round,
    HY::HydroData,
    iScen,
    t,
    NStep
    )

    path=case.DataPath
    cd(path)
    f=open("Water_volumes_levels.dat")
    line=readline(f)

    line = readline(f)
    items = split(line, " ")
    NMod = parse(Int, items[1]) #set number of modules
    water_volumes_file=zeros(Float64,HY.NMod,21);
    water_levels_file=zeros(Float64,HY.NMod,21);
    #max_head=zeros(Float64,HY.NMod);
    #min_head=zeros(Float64,HY.NMod);
    intermediate_head=zeros(Float64,HY.NMod);
    NVolumes=zeros(NMod);
  
    for iMod=1:NMod
        line = readline(f)
        items = split(line, " ")
        NVolumes[iMod] = parse(Int, items[1])
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

    for iLine = 1:2
        line = readline(f)
    end

    #=for iMod=1:NMod
        line = readline(f)                              
        max_head[iMod]=parse(Float64, strip(line))            
    end
    max_head;=# 

    #=for iMod=1:NMod
        line = readline(f)                              
        min_head[iMod]=parse(Float64, strip(line))            
    end
    min_head;=# 

    for iMod=1:NMod
        line = readline(f)                              
        intermediate_head[iMod]=parse(Float64, strip(line))            
    end
    intermediate_head; 

    #EVALUATE THE WATER LEVELS, GIVEN THE WATER VOLUMES IN THE RESERVOIR
    Level = zeros(HY.NMod)
    Head_upper = 0
    Head_lower = 0
   
    # CALCULATES THE WATER LEVELS (m a.s.l) AND THE HEAD FROM THE VOLUME RESULTS
    
    for iMod=1:HY.NMod

        for n=1:Int(NVolumes[iMod])-1
            
            if iScen == 1
                if t == 1
                    for iStep = 1:NStep
                        if iStep == 1
                            if HY.ResInit0[iMod] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n]
                            elseif HY.ResInit0[iMod] > water_volumes_file[iMod,n] && HY.ResInit0[iMod] < water_volumes_file[iMod,n+1]
                                Level[iMod] =(water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(HY.ResInit0[iMod]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if HY.ResInit0[iMod] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        else
                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n] 
                            elseif Reservoir_round[iMod,iScen,t,iStep-1]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen,t,iStep-1]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen,t,iStep-1]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        end
                            
                    end
                else
                    for iStep = 1:NStep
                        if iStep == 1        
                            if Reservoir_round[iMod,iScen,t-1,NStep] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n] 
                            elseif Reservoir_round[iMod,iScen,t-1,NStep]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen,t-1,NStep]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen,t-1,NStep]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen,t-1,NStep] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        else
                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n] 
                            elseif Reservoir_round[iMod,iScen,t,iStep-1]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen,t,iStep-1]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen,t,iStep-1]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        end
                    end
                end
            else
                if t == 1
                    for iStep = 1:NStep
                        if iStep == 1
                            if Reservoir_round[iMod,iScen-1,end,NStep] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n]
                            elseif Reservoir_round[iMod,iScen-1,end,NStep]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen-1,end,NStep]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen-1,end,NStep]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen-1,end,NStep] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        else
                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n] 
                            elseif Reservoir_round[iMod,iScen,t,iStep-1]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen,t,iStep-1]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen,t,iStep-1]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        end
                    end

                else
                    for iStep = 1:NStep
                        if iStep == 1
                            if Reservoir_round[iMod,iScen,t-1,NStep] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n] 
                            elseif Reservoir_round[iMod,iScen,t-1,NStep]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen,t-1,NStep]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen,t-1,NStep]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen,t-1,NStep] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        else
                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,n]
                                Level[iMod] = water_levels_file[iMod,n] 
                            elseif Reservoir_round[iMod,iScen,t,iStep-1]> water_volumes_file[iMod,n] && Reservoir_round[iMod,iScen,t,iStep-1]< water_volumes_file[iMod,n+1]
                                Level[iMod] = (water_levels_file[iMod,n+1]-water_levels_file[iMod,n])/(water_volumes_file[iMod,n+1]-water_volumes_file[iMod,n])*(Reservoir_round[iMod,iScen,t,iStep-1]-water_volumes_file[iMod,n])+water_levels_file[iMod,n]
                            end

                            if Reservoir_round[iMod,iScen,t,iStep-1] == water_volumes_file[iMod,Int(NVolumes[iMod])] 
                                Level[iMod] = water_levels_file[iMod,Int(NVolumes[iMod])]
                            end
                        end
                    end

                end
            end

        end
    
    end

    Head_upper = Level[1] - Level[2]
    Head_lower = Level[2] - 520
    
    return Head_data(water_volumes_file,water_levels_file,NVolumes,Head_upper,Head_lower,intermediate_head)#min_head)#max_head)

end

#SIMULATION WITH MAX HEAD CURVE

#=function efficiency_evaluation(HY::HydroData, Head::Head_data)

    @unpack (NMod,Eff,PowMaxSegTurb,DisPointTurb,PowMaxSegPump,DisPointPump,K_1_max,K_2_max,K_3_max,K_4_max,K_pump_max) = HY
    @unpack (Head_upper,Head_lower,max_head) = Head 

    P_1_1 = zeros(HY.NMod)
    P_2_1 = zeros(HY.NMod)
    K_1 = zeros(HY.NMod)
    K_2 = zeros(HY.NMod)
    K_3 = zeros(HY.NMod)
    K_4 = zeros(HY.NMod)
    Delta_Power = zeros(HY.NMod)
    

    P_1_1_pump = 0
    P_2_1_pump = 0
    K_pump = 0
    Delta_Power_pump = 0

    # Upper reservoir

    if Head_upper == max_head[1] 
        K_pump = HY.K_pump_max
    else
        eta_pump = (max_head[1] * 9810 * HY.DisPointPump[1]) / (HY.PowMaxSegPump[1] * 1000000)
        P_1_1_pump = HY.PowMaxSegPump[1]
        P_2_1_pump = (HY.DisPointPump[1] * Head_upper * 9810) / (eta_pump * 1000000)
        Delta_Power_pump = P_1_1_pump - P_2_1_pump
        K_pump = HY.K_pump_max - Delta_Power_pump
    end

    for iMod = 1:HY.NMod

        if iMod == 1

            if Head_upper == max_head[iMod] 
                K_1[iMod] = HY.K_1_max[iMod]
                K_2[iMod] = HY.K_2_max[iMod]
                K_3[iMod] = HY.K_3_max[iMod]
                K_4[iMod] = HY.K_4_max[iMod]
            else
                eta = HY.PowMaxSegTurb[iMod, 1] / (1000000 * HY.DisPointTurb[iMod, 1] * max_head[iMod] * 9810)
                P_1_1[iMod] = HY.PowMaxSegTurb[iMod, 1]
                P_2_1[iMod] = (1000000 * HY.DisPointTurb[iMod, 1] * eta * Head_upper * 9810)
                Delta_Power[iMod] = P_1_1[iMod] - P_2_1[iMod]
                K_1[iMod] = HY.K_1_max[iMod] - Delta_Power[iMod]
                K_2[iMod] = HY.K_2_max[iMod] - Delta_Power[iMod]
                K_3[iMod] = HY.K_3_max[iMod] - Delta_Power[iMod]
                K_4[iMod] = HY.K_4_max[iMod] - Delta_Power[iMod]
            end

    # Lower reservoir
    
        else
            if Head_lower == max_head[iMod] 
                K_1[iMod] = HY.K_1_max[iMod]
                K_2[iMod] = HY.K_2_max[iMod]
                K_3[iMod] = HY.K_3_max[iMod]
                K_4[iMod] = HY.K_4_max[iMod]
            else 
                eta = HY.PowMaxSegTurb[iMod, 1] / (1000000 * HY.DisPointTurb[iMod, 1] * max_head[iMod] * 9810)
                P_1_1[iMod] = HY.PowMaxSegTurb[iMod, 1]
                P_2_1[iMod] = (1000000 * HY.DisPointTurb[iMod, 1] * eta * Head_lower * 9810)
                Delta_Power[iMod] = P_1_1[iMod] - P_2_1[iMod]
                K_1[iMod] = HY.K_1_max[iMod] - Delta_Power[iMod]
                K_2[iMod] = HY.K_2_max[iMod] - Delta_Power[iMod]
                K_3[iMod] = HY.K_3_max[iMod] - Delta_Power[iMod]
                K_4[iMod] = HY.K_4_max[iMod] - Delta_Power[iMod]
            end
        end    
    end
    
    return Coeff_data(K_1, K_2, K_3, K_4, K_pump)

end =#

#SIMULATION WITH MIN HEAD CURVE

#=function efficiency_evaluation(HY::HydroData, Head::Head_data)

    @unpack (NMod,Eff,PowMaxSegTurb,DisPointTurb,PowMaxSegPump,DisPointPump,K_1_max,K_2_max,K_3_max,K_4_max,K_pump_max) = HY
    @unpack (Head_upper,Head_lower,min_head) = Head 

    P_1_1 = zeros(HY.NMod)
    P_2_1 = zeros(HY.NMod)
    K_1 = zeros(HY.NMod)
    K_2 = zeros(HY.NMod)
    K_3 = zeros(HY.NMod)
    K_4 = zeros(HY.NMod)
    Delta_Power = zeros(HY.NMod)
    

    P_1_1_pump = 0
    P_2_1_pump = 0
    K_pump = 0
    Delta_Power_pump = 0

    # Upper reservoir

    if Head_upper == min_head[1] 
        K_pump = HY.K_pump_max
    else
        eta_pump = (min_head[1] * 9810 * HY.DisPointPump[1]) / (HY.PowMaxSegPump[1] * 1000000)
        P_1_1_pump = HY.PowMaxSegPump[1]
        P_2_1_pump = (HY.DisPointPump[1] * Head_upper * 9810) / (eta_pump * 1000000)
        Delta_Power_pump = P_2_1_pump - P_1_1_pump 
        K_pump = HY.K_pump_max + Delta_Power_pump
    end

    for iMod = 1:HY.NMod

        if iMod == 1

            if Head_upper == min_head[iMod] 
                K_1[iMod] = HY.K_1_max[iMod]
                K_2[iMod] = HY.K_2_max[iMod]
                K_3[iMod] = HY.K_3_max[iMod]
                K_4[iMod] = HY.K_4_max[iMod]
            else
                eta = HY.PowMaxSegTurb[iMod, 1] / (1000000 * HY.DisPointTurb[iMod, 1] * min_head[iMod] * 9810)
                P_1_1[iMod] = HY.PowMaxSegTurb[iMod, 1]
                P_2_1[iMod] = (1000000 * HY.DisPointTurb[iMod, 1] * eta * Head_upper * 9810)
                Delta_Power[iMod] = P_2_1[iMod] - P_1_1[iMod] 
                K_1[iMod] = HY.K_1_max[iMod] + Delta_Power[iMod]
                K_2[iMod] = HY.K_2_max[iMod] + Delta_Power[iMod]
                K_3[iMod] = HY.K_3_max[iMod] + Delta_Power[iMod]
                K_4[iMod] = HY.K_4_max[iMod] + Delta_Power[iMod]
            end

    # Lower reservoir
    
        else
            if Head_lower == min_head[iMod] 
                K_1[iMod] = HY.K_1_max[iMod]
                K_2[iMod] = HY.K_2_max[iMod]
                K_3[iMod] = HY.K_3_max[iMod]
                K_4[iMod] = HY.K_4_max[iMod]
            else 
                eta = HY.PowMaxSegTurb[iMod, 1] / (1000000 * HY.DisPointTurb[iMod, 1] * min_head[iMod] * 9810)
                P_1_1[iMod] = HY.PowMaxSegTurb[iMod, 1]
                P_2_1[iMod] = (1000000 * HY.DisPointTurb[iMod, 1] * eta * Head_lower * 9810)
                Delta_Power[iMod] = P_2_1[iMod] - P_1_1[iMod]  
                K_1[iMod] = HY.K_1_max[iMod] + Delta_Power[iMod]
                K_2[iMod] = HY.K_2_max[iMod] + Delta_Power[iMod]
                K_3[iMod] = HY.K_3_max[iMod] + Delta_Power[iMod]
                K_4[iMod] = HY.K_4_max[iMod] + Delta_Power[iMod]
            end
        end    
    end
    
    return Coeff_data(K_1, K_2, K_3, K_4, K_pump)

end=# 

#SIMULATION WITH INTERMEDIATE HEAD CURVE

function efficiency_evaluation(HY::HydroData, Head::Head_data)

    @unpack (NMod,Eff,PowMaxSegTurb,DisPointTurb,PowMaxSegPump,DisPointPump,K_1_max,K_2_max,K_3_max,K_4_max,K_pump_max) = HY
    @unpack (Head_upper,Head_lower,intermediate_head) = Head 

    P_1_1 = zeros(HY.NMod)
    P_2_1 = zeros(HY.NMod)
    K_1 = zeros(HY.NMod)
    K_2 = zeros(HY.NMod)
    K_3 = zeros(HY.NMod)
    K_4 = zeros(HY.NMod)
    Delta_Power = zeros(HY.NMod)
    

    P_1_1_pump = 0
    P_2_1_pump = 0
    K_pump = 0
    Delta_Power_pump = 0

    # Upper reservoir

    if Head_upper == intermediate_head[1] 
        K_pump = HY.K_pump_max
    else
        eta_pump = (intermediate_head[1] * 9810 * HY.DisPointPump[1]) / (HY.PowMaxSegPump[1] * 1000000)
        P_1_1_pump = HY.PowMaxSegPump[1]
        P_2_1_pump = (HY.DisPointPump[1] * Head_upper * 9810) / (eta_pump * 1000000)
        if P_1_1_pump > P_2_1_pump
            Delta_Power_pump = P_1_1_pump - P_2_1_pump 
            K_pump = HY.K_pump_max - Delta_Power_pump
        elseif P_1_1_pump < P_2_1_pump
            Delta_Power_pump = P_2_1_pump - P_1_1_pump 
            K_pump = HY.K_pump_max + Delta_Power_pump
        end 
    end

    for iMod = 1:HY.NMod

        if iMod == 1

            if Head_upper == intermediate_head[iMod] 
                K_1[iMod] = HY.K_1_max[iMod]
                K_2[iMod] = HY.K_2_max[iMod]
                K_3[iMod] = HY.K_3_max[iMod]
                K_4[iMod] = HY.K_4_max[iMod]
            else
                eta = HY.PowMaxSegTurb[iMod, 1] / (1000000 * HY.DisPointTurb[iMod, 1] * intermediate_head[iMod] * 9810)
                P_1_1[iMod] = HY.PowMaxSegTurb[iMod, 1]
                P_2_1[iMod] = (1000000 * HY.DisPointTurb[iMod, 1] * eta * Head_upper * 9810)
                if P_1_1_pump > P_2_1_pump
                    Delta_Power[iMod] = P_1_1[iMod] - P_2_1[iMod]
                    K_1[iMod] = HY.K_1_max[iMod] - Delta_Power[iMod]
                    K_2[iMod] = HY.K_2_max[iMod] - Delta_Power[iMod]
                    K_3[iMod] = HY.K_3_max[iMod] - Delta_Power[iMod]
                    K_4[iMod] = HY.K_4_max[iMod] - Delta_Power[iMod]
                elseif P_1_1_pump < P_2_1_pump
                    Delta_Power[iMod] = P_2_1[iMod] - P_1_1[iMod] 
                    K_1[iMod] = HY.K_1_max[iMod] + Delta_Power[iMod]
                    K_2[iMod] = HY.K_2_max[iMod] + Delta_Power[iMod]
                    K_3[iMod] = HY.K_3_max[iMod] + Delta_Power[iMod]
                    K_4[iMod] = HY.K_4_max[iMod] + Delta_Power[iMod]
                end
            end

    # Lower reservoir
    
        else
            if Head_lower == intermediate_head[iMod] 
                K_1[iMod] = HY.K_1_max[iMod]
                K_2[iMod] = HY.K_2_max[iMod]
                K_3[iMod] = HY.K_3_max[iMod]
                K_4[iMod] = HY.K_4_max[iMod]
            else 
                eta = HY.PowMaxSegTurb[iMod, 1] / (1000000 * HY.DisPointTurb[iMod, 1] * intermediate_head[iMod] * 9810)
                P_1_1[iMod] = HY.PowMaxSegTurb[iMod, 1]
                P_2_1[iMod] = (1000000 * HY.DisPointTurb[iMod, 1] * eta * Head_lower * 9810)
                if P_1_1_pump > P_2_1_pump
                    Delta_Power[iMod] = P_1_1[iMod] - P_2_1[iMod]
                    K_1[iMod] = HY.K_1_max[iMod] - Delta_Power[iMod]
                    K_2[iMod] = HY.K_2_max[iMod] - Delta_Power[iMod]
                    K_3[iMod] = HY.K_3_max[iMod] - Delta_Power[iMod]
                    K_4[iMod] = HY.K_4_max[iMod] - Delta_Power[iMod]
                elseif P_1_1_pump < P_2_1_pump
                    Delta_Power[iMod] = P_2_1[iMod] - P_1_1[iMod] 
                    K_1[iMod] = HY.K_1_max[iMod] + Delta_Power[iMod]
                    K_2[iMod] = HY.K_2_max[iMod] + Delta_Power[iMod]
                    K_3[iMod] = HY.K_3_max[iMod] + Delta_Power[iMod]
                    K_4[iMod] = HY.K_4_max[iMod] + Delta_Power[iMod]
                end
            end
        end    
    end
    
    return Coeff_data(K_1, K_2, K_3, K_4, K_pump)

end 
