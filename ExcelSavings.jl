
using DataFrames,XLSX

function data_saving(runMode, SimScen, InputParameters)                                                                                                                                               # Genero la cartella

    hour=string(now())
    a=replace(hour,':'=>'-')

    @unpack (NStep, NSeg, NSimScen, NStage, NHoursStep ) = InputParameters
    @unpack (scenarios, scenStates) = SimScen

    NSeg_up = NSeg[1]
    NSeg_low = NSeg[2]

    type=""
    if HY.Station_with_pump==0
        type="Trad"
    else
        type="Pump"
    end

    envConst=""
    if runMode.envConst
        envConst="_wEnv"
    else
        envConst="_woEnv"
    end

    rampConst=""
    if runMode.ramping_constraints
        rampConst="_wRamp"
    else
        rampConst="_woRamp"
    end

    mode="$type,$envConst,$rampConst, $NStage stages, $NStep steps, $NSeg_up upper, $NSeg_low lower, $a"

    folder2 = "$mode"                                                               # Sto creando una nuova cartella "Scenario_" all'interno della cartella RESULTS
    mkdir(folder2)                                                                                                                                                             # Genero la cartella
    cd(folder2)
    main=pwd()
    println("Folder in which data are saved:", main)

    # SAVING DATA FOR WEEK 

    folder= "Weeks"
    mkdir(folder)
    cd(folder)

    for iWeek=1:NStage

        inflows_upper = DataFrame()
        inflows_lower = DataFrame()
        prices = DataFrame()
        states = DataFrame()
        ResInit= DataFrame()
        statesSim = zeros(NSimScen)
        prod_Uturbine = DataFrame()
        prod_Lturbine = DataFrame()
        pump = DataFrame()
        revenue_Uturbine = DataFrame()
        revenue_Lturbine = DataFrame()
        cost_pump = DataFrame()
        discharge_UpTurbine = DataFrame()
        discharge_loTurbine = DataFrame()
        discharge_pump = DataFrame()
        upper_volume = DataFrame()
        lower_volume = DataFrame()
        Level_variations_upper= DataFrame()         
        Level_variations_lower = DataFrame()
        Spillage_upper = DataFrame()
        Spillage_lower = DataFrame()

        prod_turbine = zeros(HY.NMod,NSimScen,NStage,NStep)
        prod_pump = zeros(HY.NMod,NSimScen,NStage,NStep)

        for iMod=1:HY.NMod
            for iScen= 1:NSimScen
                for iStage=1:NStage
                    for iStep=1:NStep
                        prod_turbine[iMod,iScen,iStage,iStep] = ResultsSim.Production[iMod,iScen,iStage,iStep]*NHoursStep   #MWh
                        prod_pump[iMod,iScen,iStage,iStep] = ResultsSim.Pumping[iMod,iScen,iStage,iStep]*NHoursStep
                    end
                end
            end
        end

        for iScen=1:NSimScen
                inflows_upper[!,"Scen_$iScen"] = ResultsSim.inflow[1,iScen,iWeek,:]
                inflows_lower[!,"Scen_$iScen"] = ResultsSim.inflow[2,iScen,iWeek,:]
                prices[!,"Scen_$iScen"] = ResultsSim.price[1,iScen,iWeek,:]
                statesSim[iScen] = SimScen.scenStates[iScen][iWeek]
                prod_Uturbine[!,"Scen_$iScen"] = prod_turbine[1,iScen,iWeek,:]
                prod_Lturbine[!,"Scen_$iScen"] = prod_turbine[2,iScen,iWeek,:]
                pump[!,"Scen_$iScen"] = prod_pump[1,iScen,iWeek,:]
                revenue_Uturbine[!,"Scen_$iScen"] = ResultsSim.turbine_profit_timestep[1,iScen,iWeek,:]
                revenue_Lturbine[!,"Scen_$iScen"] = ResultsSim.turbine_profit_timestep[2,iScen,iWeek,:]
                cost_pump[!,"Scen_$iScen"] = ResultsSim.pumping_costs_timestep[1,iScen,iWeek,:]
                discharge_UpTurbine[!,"Scen_$iScen"] = ResultsSim.totDischarge[1,iScen,iWeek,:]
                discharge_loTurbine[!,"Scen_$iScen"] = ResultsSim.totDischarge[2,iScen,iWeek,:]
                discharge_pump[!,"Scen_$iScen"] = ResultsSim.totPumped[iScen,iWeek,:]
                upper_volume[!,"Scen_$iScen"] = ResultsSim.Reservoir[1,iScen,iWeek,:]
                lower_volume[!,"Scen_$iScen"] = ResultsSim.Reservoir[2,iScen,iWeek,:]
                Level_variations_upper[!,"Scen_$iScen"] = Results_Water_levels.water_level_variations[1,iScen,iWeek,:]
                Level_variations_lower[!,"Scen_$iScen"] = Results_Water_levels.water_level_variations[2,iScen,iWeek,:]
                Spillage_upper[!,"Scen_$iScen"] = ResultsSim.Spillage[1,iScen,iWeek,:]
                Spillage_lower[!,"Scen_$iScen"] = ResultsSim.Spillage[2,iScen,iWeek,:]

            #=if iWeek==1      
                if iScen==1 
                    ResInit[!,"Scen_$iScen"] = HY.ResInit0[:]
                else    #iScen>1
                    ResInit[!,"Scen_$iScen"] = ResultsSim.Reservoir[:,iScen-1,end,end]
                end
            else
                ResInit[!,"Scen_$iScen"] = ResultsSim.Reservoir[:,iScen,iWeek-1,end]
            end =#

            if iScen==1      
                if iWeek==1 
                    ResInit[!,"Scen_$iScen"] = HY.ResInit0[:]
                else    #iWeek>1
                    ResInit[!,"Scen_$iScen"] = ResultsSim.Reservoir[:,iScen,iWeek-1,end]
                end
            else # iScen>1
                if iWeek==1 
                    ResInit[!,"Scen_$iScen"] = ResultsSim.Reservoir[:,iScen-1,end,end]
                else    #iWeek>1
                    ResInit[!,"Scen_$iScen"] = ResultsSim.Reservoir[:,iScen,iWeek-1,end]
                end
            end
            
        end
        
        states[!,"Week_$iWeek"] = statesSim[:]

        XLSX.writetable("Week_$iWeek.xlsx",overwrite=true,
        ResInit_Mm3=(collect(DataFrames.eachcol(ResInit)),DataFrames.names(ResInit)),
        State=(collect(DataFrames.eachcol(states)),DataFrames.names(states)),
        Inflows_upper_Mm3=(collect(DataFrames.eachcol(inflows_upper)),DataFrames.names(inflows_upper)),
        Inflows_lower_Mm3=(collect(DataFrames.eachcol(inflows_lower)),DataFrames.names(inflows_lower)),
        Prices_€alMWh=(collect(DataFrames.eachcol(prices)),DataFrames.names(prices)),
        Prod_Uturbine_MWh = (collect(DataFrames.eachcol(prod_Uturbine)),DataFrames.names(prod_Uturbine)),
        Prod_Lturbine_MWh = (collect(DataFrames.eachcol(prod_Lturbine)),DataFrames.names(prod_Lturbine)),
        Pumping_MWh = (collect(DataFrames.eachcol(pump)),DataFrames.names(pump)),
        Rev_Uturbine_€ = (collect(DataFrames.eachcol(revenue_Uturbine)),DataFrames.names(revenue_Uturbine)),
        Rev_Lturbine_€ = (collect(DataFrames.eachcol(revenue_Lturbine)),DataFrames.names(revenue_Lturbine)),
        Cost_pump_€ = (collect(DataFrames.eachcol(cost_pump)),DataFrames.names(cost_pump)),
        dis_UpTurbine_m3s = (collect(DataFrames.eachcol(discharge_UpTurbine)),DataFrames.names(discharge_UpTurbine)),
        dis_LoTurbine_m3s = (collect(DataFrames.eachcol(discharge_loTurbine)),DataFrames.names(discharge_loTurbine)),
        dis_pump_m3s = (collect(DataFrames.eachcol(discharge_pump)),DataFrames.names(discharge_pump)),
        upper_volume_Mm3 = (collect(DataFrames.eachcol(upper_volume)),DataFrames.names(upper_volume)),
        lower_volume_Mm3 = (collect(DataFrames.eachcol(lower_volume)),DataFrames.names(lower_volume)),
        water_var_Upper_m = (collect(DataFrames.eachcol(Level_variations_upper)),DataFrames.names(Level_variations_upper)),
        water_var_Lower_m =(collect(DataFrames.eachcol(Level_variations_lower)),DataFrames.names(Level_variations_lower)),
        Spillage_upper_m3s =(collect(DataFrames.eachcol(Spillage_upper)),DataFrames.names(Spillage_upper)),
        Spillage_lower_m3s =(collect(DataFrames.eachcol(Spillage_lower)),DataFrames.names(Spillage_lower))
        )

    end 

    cd(main)
   
    #SAVING PRODUCTION AND REVEUE VALUES

    prod= DataFrame();
    Production_factor_upper = DataFrame()
    Production_factor_lower = DataFrame()
    Production_factor_pump = DataFrame()

    total_turbine_production=zeros(HY.NMod,NSimScen)
    total_pump_power = zeros(HY.NMod,NSimScen)
    total_profit_turbine = zeros(HY.NMod,NSimScen)
    total_cost_pump = zeros(HY.NMod,NSimScen)
        
    for iScen=1:NSimScen
        for iMod=1:HY.NMod
            for iStage=1:NStage
                for iStep=1:NStep
                    total_turbine_production[iMod,iScen]=total_turbine_production[iMod,iScen]+ResultsSim.Production[iMod,iScen,iStage,iStep]*NHoursStep
                    total_pump_power[iMod,iScen] = total_pump_power[iMod,iScen]+ ResultsSim.Pumping[iMod,iScen,iStage,iStep]*NHoursStep
                    total_profit_turbine[iMod,iScen] = total_profit_turbine[iMod,iScen]+ ResultsSim.turbine_profit_timestep[iMod,iScen,iStage,iStep]
                    total_cost_pump[iMod,iScen]= total_cost_pump[iMod,iScen]+ResultsSim.pumping_costs_timestep[iMod,iScen,iStage,iStep]
                end    
            end
        end
    end
        
    prod[!,"Scenario"]= 1:1:NSimScen;

    prod[!,"Production_turbine_upper [MWh]"]= total_turbine_production[1,:];
    prod[!,"Profit_turbine_upper [€]"] = total_profit_turbine[1,:];
    prod[!,"Power_for_pump [MWh]"] = total_pump_power[1,:];
    prod[!,"Cost for pump [€]"] = total_cost_pump[1,:];
    prod[!,"Production_turbine_lower [MWh]"] = total_turbine_production[2,:];
    prod[!,"Profit_turbine_lower [€]"] = total_profit_turbine[2,:];
    

    Production_factor_upper[!,"Scenario"]=1:1:NSimScen;
    Production_factor_lower[!,"Scenario"]=1:1:NSimScen;
    Production_factor_pump[!,"Scenario"]=1:1:NSimScen;

    range1=["pf==0","0<pf<0.25","0.25<=pf<=0.5","0.5<pf<=0.75","0.75<pf<1","pf==1","pf>1"]

    for i=1:length(range1)
        Production_factor_upper[!,range1[i]] = Production_factors.nsteps_turbines[1,:,i]
        Production_factor_lower[!,range1[i]]= Production_factors.nsteps_turbines[2,:,i]
        Production_factor_pump[!,range1[i]] = Production_factors.nsteps_pumps[1,:,i]
    end

    XLSX.writetable("Yearly_production.xlsx",overwrite=true,
    Production=(collect(DataFrames.eachcol(prod)),DataFrames.names(prod)),
    Production_factor_upper=(collect(DataFrames.eachcol(Production_factor_upper)),DataFrames.names(Production_factor_upper)),
    Production_factor_lower=(collect(DataFrames.eachcol(Production_factor_lower)),DataFrames.names(Production_factor_lower)),
    Production_factor_pump=(collect(DataFrames.eachcol(Production_factor_pump)),DataFrames.names(Production_factor_pump))
    )

    cd(main)

    #SAVING WATER LEVEL VARIATIONS FREQUENCIES

    frequency_allScenarios_upper=DataFrame()
    frequency_allScenarios_lower=DataFrame()
    frequency_allScenarios_upper[!,"Scenarios"] = 1:1:NSimScen
    frequency_allScenarios_lower[!,"Scenarios"] = 1:1:NSimScen
    
    ranges=["v<-0.20","-0.20<=v<-0.15","-0.15<=v<-0.10","-0.10<v<-0.05","-0.05<=v<0.00","v==0","0<v<=0.05","0.05<v<=0.10","0.10<v<=0.15","0.15<v<=0.20","v>0.20"]

    for i=1:length(ranges)
        frequency_allScenarios_upper[!,ranges[i]] = Results_Water_levels.frequency[1,:,i]
        frequency_allScenarios_lower[!,ranges[i]] = Results_Water_levels.frequency[2,:,i]
    end

    XLSX.writetable("Frequency_Water_levels.xlsx",overwrite=true,
    Upper_level_frequencies=(collect(DataFrames.eachcol(frequency_allScenarios_upper)),DataFrames.names(frequency_allScenarios_upper)),
    Lower_level_frequencies=(collect(DataFrames.eachcol(frequency_allScenarios_lower)),DataFrames.names(frequency_allScenarios_lower)),
    )
 
    cd(main)        
    
    #SAVING WATER VALUES

#=    WaterValuesFolder= "Water Values"
    mkdir(WaterValuesFolder)
    cd(WaterValuesFolder)

    for iSeg=1:length(ResSeg)

        Water_Values_Upper = DataFrame()
        Water_Values_Lower = DataFrame()
        Alpha_Table = DataFrame()

        for iState=1:NStates
            Water_Values_Upper[!,"State_$iState"]= ResultsSDP.WVTable[:,iState,iSeg,1]
            Water_Values_Lower[!,"State_$iState"]= ResultsSDP.WVTable[:,iState,iSeg,2]
            Alpha_Table[!,"State_$iState"] = ResultsSDP.AlphaTable[:,iState,iSeg]
        end

        Upper = ResultsSDP.ResSeg[iSeg][1]
        lower = ResultsSDP.ResSeg[iSeg][2]

        XLSX.writetable("C$iSeg - Upper $Upper - Lower $lower .xlsx",overwrite=true,
            Water_values_upper_€_Mm3 = (collect(DataFrames.eachcol(Water_Values_Upper)),DataFrames.names(Water_Values_Upper)),
            Water_values_lower_€_Mm3 = (collect(DataFrames.eachcol(Water_Values_Lower)),DataFrames.names(Water_Values_Lower)),
            Alpha_values_€ = (collect(DataFrames.eachcol(Alpha_Table)),DataFrames.names(Alpha_Table))
        )
    end 
     
    cd(main)

    # SAVING ERROR EVALUATION

    if runMode.error_evaluation 
  
        path="Error"
        mkdir(path)
        cd(path)
      
        error_evaluation(InputParameters)
      
    end=#

    return ("Saved data")
end





