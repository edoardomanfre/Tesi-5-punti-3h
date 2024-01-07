

positive_var=0
negative_var=0
min_flow_var=0

for iScen=1:NSimScen
    for iMod=1:HY.NMod
        for iStage=1:NStage
            for iStep =1:InputParameters.NStep

                if ResultsSim.Res_slack_pos[iMod,iScen,iStage,iStep] > 0
                    println("Slack_positive_variation_Module:",iMod," Scenario:",iScen," Stage:",iStage," iStep:",iStep," Value:",ResultsSim.Res_slack_pos[iMod,iScen,iStage,iStep])
                    positive_var = positive_var+1
                end

            end
        end
    end
end

println("N.times_Slack_Positive_var:",positive_var)
println("N.times_Slack_Negative_var:",negative_var)
println("N.times_Slack_Min_env:",min_flow_var)

negative_var=0

for iScen=1:NSimScen
   # for iMod=1:HY.NMod
        for iStage=1:NStage
            for iStep =1:InputParameters.NStep

                if Results_Water_levels.water_level_variations[1,iScen,iStage,iStep] < - 0.55
                    println("Negative variations: "," Scenario:",iScen," Stage:",iStage," iStep:",iStep, " Variation:",Results_Water_levels.water_level_variations[1,iScen,iStage,iStep])
                    negative_var = negative_var+1
                end

            end
        end
   # end
end
negative_var

pos=0
for iScen=1:NSimScen
   # for iMod=1:HY.NMod
        for iStage=1:22
            for iStep =1:InputParameters.NStep

                if Results_Water_levels.water_level_variations[1,iScen,iStage,iStep] > 0.55
                    println("Positive variations: 1"," Scenario:",iScen," Stage:",iStage," iStep:",iStep, " Variation:",Results_Water_levels.water_level_variations[1,iScen,iStage,iStep])
                    println("Reservoir volueme: ", ResultsSim.Reservoir[1,iScen,iStage,iStep])
                    println("Positive slack variable: ", ResultsSim.Res_slack_pos[1,iScen,iStage,iStep])
                    println("Inflow: ", ResultsSim.inflow[1,iScen,iStage,iStep])
                    pos = pos+1
                end

            end
        end
   # end
end
pos



#for iMod=1:HY.NMod
# SPILLAGE FROM UPPER RESERVOIR
spill=0
    for iScen=1:NSimScen
        for iStage=1:NStage
            for iStep =1:InputParameters.NStep

                if ResultsSim.Spillage[1,iScen,iStage,iStep] > 0
                    println("Spillage from:1 ", "Scenario:",iScen," Stage:",iStage," iStep:",iStep, " Water spilled:",ResultsSim.Spillage[1,iScen,iStage,iStep])
                    println("Inflow for:1 "," Scenario:",iScen," Stage:", iStage," iStep",iStep," Inflow: ",ResultsSim.inflow[1,iScen,iStage,iStep])
                    spill = spill+1
                end
            end   
        end
    end
#end
println("Times Spillage from upper, $spill")

# SPILLAGE FROM LOWER
spill2=0
for iScen=1:NSimScen
    for iStage=1:NStage
        for iStep =1:InputParameters.NStep

            if ResultsSim.Spillage[2,iScen,iStage,iStep] > 0
                println("Spillage from:2 ", "Scenario:",iScen," Stage:",iStage," iStep:",iStep, " Water spilled:",ResultsSim.Spillage[2,iScen,iStage,iStep])
                println("Inflow for:2 "," Scenario:",iScen," Stage:", iStage," iStep",iStep," Inflow: ",ResultsSim.inflow[2,iScen,iStage,iStep])
                spill2 = spill2+1
            end
        end   
    end
end
#end
println("Times Spillage from upper, $spill")

using DataFrames,XLSX

    frequency_allScenarios_upper=DataFrame()
    frequency_allScenarios_lower=DataFrame()
    frequency_allScenarios_upper[!,"Scenarios"] = 1:1:NSimScen
    frequency_allScenarios_lower[!,"Scenarios"] = 1:1:NSimScen
    
    ranges=["v<-0.50","-0.5<=v<-0.30","-0.30<=v<-0.10","-0.10<v<-0.05","-0.05<=v<0.00","v==0","0<v<=0.05","0.05<v<=0.10","0.1<v<=0.30","0.30<v<=0.5","v>0.5"]

    for i=1:length(ranges)
        frequency_allScenarios_upper[!,ranges[i]] = Results_Water_levels.frequency[1,:,i]
        frequency_allScenarios_lower[!,ranges[i]] = Results_Water_levels.frequency[2,:,i]
    end

   
    XLSX.writetable("3_Water Frequency.xlsx",overwrite=true,
    Upper_level_frequencies=(collect(DataFrames.eachcol(frequency_allScenarios_upper)),DataFrames.names(frequency_allScenarios_upper)),
    Lower_level_frequencies=(collect(DataFrames.eachcol(frequency_allScenarios_lower)),DataFrames.names(frequency_allScenarios_lower)),
    )


    # FIND MAXIMUM SPILLED WATER IN THAT week
    for iScen=1:NSimScen
        for iStage=1:NStage
           # for iMod=1:HY.NMod
                for iStep=1:7
                    if ResultsSim.Spillage[1,iScen,iStage,iStep] >0
                        a= findmax(ResultsSim.Spillage[1,iScen,iStage,:])
                        println("Maximum spiallage for module 1, scenario $iScen, week $iStage " ,a)
                        
                        b=findmax(ResultsSim.inflow[1,iScen,iStage,:])
                        println("Maximum inflow for module 1, scenario $iScen, week $iStage ",b)

                    end
                end
           # end
        end
    end


    for iScen=1:NSimScen
        for iStage=1:NStage
           # for iMod=1:HY.NMod
                for iStep=1:7
                    if ResultsSim.Spillage[2,iScen,iStage,iStep] >0
                    a= findmax(ResultsSim.Spillage[2,iScen,iStage,:])
                    println("Maximum spillage for module 2, scenario $iScen, week $iStage " ,a)
                    end
                end
           # end
        end
    end



    production_factors_UPPER_turbines=DataFrame();
    production_factor_LOWER_turbines = DataFrame();
    production_factor_pump = DataFrame();
    net_production= DataFrame();
    mean_days=DataFrame();
    total_net_production=zeros(HY.NMod,NSimScen);
    
    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            total_net_production[iMod,iScen]=sum(ResultsSim.Net_production[iMod,iScen,:,:])
        end
    end
    
    net_production[!,"Scenario"]= 1:1:NSimScen;
    net_production[!,"Net_production_upper [MW]"]= total_net_production[1,:];
    net_production[!,"Net_profit_upper [€]"] = ResultsSim.annual_profit_each_reservoir[1,:];
    net_production[!,"Net_production_lower [MW]"] = total_net_production[2,:];
    net_production[!,"Net_profit_lower [€]"] = ResultsSim.annual_profit_each_reservoir[2,:];
  

    production_factors_UPPER_turbines[!,"Scenario"]=1:1:NSimScen;
    production_factor_LOWER_turbines[!,"Scenario"]=1:1:NSimScen;
    production_factor_pump[!,"Scenario"]=1:1:NSimScen;

    range1=["pf==0","0<pf<0.25","0.25<=pf<=0.5","0.5<pf<=0.75","0.75<pf<1","pf==1","pf>1"]

    for i=1:length(range1)
        production_factors_UPPER_turbines[!,range1[i]]=Production_factors.nsteps_turbines[1,:,i]
        production_factor_LOWER_turbines[!,range1[i]]= Production_factors.nsteps_turbines[2,:,i]
        production_factor_pump[!,range1[i]] = Production_factors.nsteps_pumps[1,:,i]
    end
  
    mean_days[!,"Production Factors"]= range1[:];
    mean_days[!,"Mean Upper Turbines"]=Production_factors.mean_nsteps_turbines[1,:];
    mean_days[!,"Mean Lower Tubines"]=Production_factors.mean_nsteps_turbines[2,:];
    mean_days[!,"Mean Pump"]= Production_factors.mean_nsteps_pumps[1,:]
  
  
    XLSX.writetable("Prod_factors.xlsx",overwrite=true,
    Production=(collect(DataFrames.eachcol(net_production)),DataFrames.names(net_production)),
    Production_factors_upper=(collect(DataFrames.eachcol(production_factors_UPPER_turbines)),DataFrames.names(production_factors_UPPER_turbines)),
    Production_factor_lower=(collect(DataFrames.eachcol(production_factor_LOWER_turbines)),DataFrames.names(production_factor_LOWER_turbines)),
    Production_factor_pump=(collect(DataFrames.eachcol(production_factor_pump)),DataFrames.names(production_factor_pump)),
    Mean_factors=(collect(DataFrames.eachcol(mean_days)),DataFrames.names(mean_days))
    )