# OUTPUT VALUES FOR GIVEN SCENARIO N.RANDOM
using Statistics
using Plots
pyplot()
using DataFrames
using XLSX
using DelimitedFiles


function SavePlots(InputParameters::InputParam)

    Scenarios=[36 53];
    scenario_max_inflow=36;
    scenario_min_inflow=53;
    println("Plots for scenarios $scenario_max_inflow, and $scenario_min_inflow")

    # HOW TO CONCATENATE MORE VECTORS - HAVE VALUES OF 365 DAYS
    NStep= InputParameters.NStep;
    concatenation_production_turbines=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_production_pump=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_discharge_turbines=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_discharge_pumps=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_inflows=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_reservoir=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_price=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_profit_turbines=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_cost_pump=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_waterLevel_variations = zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_waterLevel = zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_PF_Turbines = zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_PF_Pump = zeros(HY.NMod,NSimScen,NStep*NStage);

    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            for iStage=1:NStage
                concatenation_production_turbines[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.Production[iMod,iScen,iStage,:]*InputParameters.NHoursStep
                concatenation_production_pump[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.Pumping[iMod,iScen,iStage,:]*InputParameters.NHoursStep
                concatenation_discharge_turbines[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.totDischarge[iMod,iScen,iStage,:]
                concatenation_discharge_pumps[1,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.totPumped[iScen,iStage,:]
                concatenation_inflows[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.inflow[iMod,iScen,iStage,:]
                concatenation_reservoir[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.Reservoir[iMod,iScen,iStage,:]
                concatenation_price[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.price[iMod,iScen,iStage,:]
                concatenation_profit_turbines[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)] =ResultsSim.turbine_profit_timestep[iMod,iScen,iStage,:]*InputParameters.NHoursStep
                concatenation_cost_pump[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)] = ResultsSim.pumping_costs_timestep[iMod,iScen,iStage,:]*InputParameters.NHoursStep
                concatenation_waterLevel_variations[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)] = Results_Water_levels.water_level_variations[iMod,iScen,iStage,:]
                concatenation_waterLevel[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)] = Results_Water_levels.Water_levels_Simulation[iMod,iScen,iStage,:]
                concatenation_PF_Turbines[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)] = Production_factors.step_production_factor_turbines[iMod,iScen,iStage,:]
                concatenation_PF_Pump[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)] = Production_factors.step_factor_pump[iMod,iScen,iStage,:]
            end
        end
    end

    # SAVE DATA WET YEAR and DRY YEAR - CONCATENATED VALUES
    
    for i in Scenarios
    
        Prod_Turbines = DataFrame();
        Prod_Pump = DataFrame();
        Profit_turbine = DataFrame();
        Cost_Pump = DataFrame();
        Inflow = DataFrame();
        Reservoir_volume= DataFrame();
        Reservoir_level = DataFrame();
        Variations_water= DataFrame();
        PF_Turbine = DataFrame();
        PF_Pump = DataFrame();

        for iStep=NStage*NStep
            for iMod=1:HY.NMod
                Prod_Turbines[!,"Module_$iMod"] = concatenation_production_turbines[iMod,i,:]
                Prod_Pump[!,"Module_$iMod"]  = concatenation_production_pump[iMod,i,:]
                Profit_turbine[!,"Module_$iMod"]  = concatenation_profit_turbines[iMod,i,:]
                Cost_Pump[!,"Module_$iMod"] = concatenation_cost_pump[iMod,i,:]
                Inflow[!,"Module_$iMod"] = concatenation_inflows[iMod,i,:]
                Reservoir_volume[!,"Module_$iMod"]  = concatenation_reservoir[iMod,i,:]
                Reservoir_level[!,"Module_$iMod"]  = concatenation_waterLevel[iMod,i,:]
                Variations_water[!,"Module_$iMod"]  = concatenation_waterLevel_variations[iMod,i,:]
                PF_Turbine[!,"Module_$iMod"] = concatenation_PF_Turbines[iMod,i,:]
                PF_Pump[!,"Module_$iMod"] = concatenation_PF_Pump[iMod,i,:]
            end
        end
    
        XLSX.writetable("Scenario $i.xlsx",overwrite=true,
        Prod_Turbines_MWh = (collect(DataFrames.eachcol(Prod_Turbines)),DataFrames.names(Prod_Turbines)),
        Prod_Pump_MWh = (collect(DataFrames.eachcol(Prod_Pump)),DataFrames.names(Prod_Pump)),
        Profit_turbine€ = (collect(DataFrames.eachcol(Profit_turbine)),DataFrames.names(Profit_turbine)),
        Cost_pump€ = (collect(DataFrames.eachcol(Cost_Pump)),DataFrames.names(Cost_Pump)),
        Reservoir_volume = (collect(DataFrames.eachcol(Reservoir_volume)),DataFrames.names(Reservoir_volume)),
        Reservoir_level = (collect(DataFrames.eachcol(Reservoir_level)),DataFrames.names(Reservoir_level)),
        Variations_water = (collect(DataFrames.eachcol(Variations_water)),DataFrames.names(Variations_water)),
        PF_turbine = (collect(DataFrames.eachcol(PF_Turbine)),DataFrames.names(PF_Turbine)),
        PF_Pump = (collect(DataFrames.eachcol(PF_Pump)),DataFrames.names(PF_Pump)),
        Inflow = (collect(DataFrames.eachcol(Inflow)),DataFrames.names(Inflow))
    )

end

    x=1:1:NStep*NStage;

    #Come cambiare le scale??   MAXIMUM INFLOW
    A=plot(x,concatenation_inflows[1,scenario_max_inflow,:],size=(1000,400), c=:blue, label="Inflows upper")
    plot!(twinx(),concatenation_reservoir[1,scenario_max_inflow,:],c=:red,label="Reservoir volume upper")

    B=plot(x,concatenation_inflows[2,scenario_max_inflow,:],size=(1000,400), label="Inflows lower")
    plot!(twinx(),concatenation_reservoir[2,scenario_max_inflow,:],c=:red,label="Reservoir volume lower")

    C=plot(A,B,layout=(2,1),size=(1200,1000),title ="Scenario with maximum inflow")

    # MINIMUM INFLOW
    a=plot(x,concatenation_inflows[1,scenario_min_inflow,:],size=(1000,500), c=:blue, label="Inflows upper")
    plot!(twinx(),concatenation_reservoir[1,scenario_min_inflow,:],c=:red,label="Reservoir volume upper")

    b=plot(x,concatenation_inflows[2,scenario_min_inflow,:],size=(1000,500), label="Inflows lower")
    plot!(twinx(),concatenation_reservoir[2,scenario_min_inflow,:],c=:red,label="Reservoir volume lower")

    c=plot(a,b,layout=(2,1),size=(1200,1000),title ="Scenario with minimum inflow")    
    
    # MAXIMUM INFLOW SCENARIOS : RESERVOIRS' VOLUMES
    D=plot(x,concatenation_reservoir[1,scenario_max_inflow,:],size=(1200,700),xlabel="N steps",ylabel="Reservoir's volume [Mm3]",label="Upper reservoir")
    plot!(D,x,concatenation_reservoir[2,scenario_max_inflow,:], title ="Scenario with maximum inflow", label = "Lower reservoir")
    D

    # Water Level Variations
    F=plot(x,concatenation_waterLevel_variations[1,scenario_max_inflow,:],size=(1200,700),xlabel="N steps",ylabel="Reservoir's water level variations [Mm3]",title="Scenario with maximum inflow",label="Upper reservoir")
    plot!(F,x,concatenation_waterLevel_variations[2,scenario_max_inflow,:],label="Lower reservoir")
    F

    # Production
    T=plot(x,concatenation_production_turbines[1,scenario_max_inflow,:],size=(1200,700),xlabel="N Steps", ylabel="Production [MW]", title="Scenario with maximum inflow",label="Upper reservoir")
    plot!(T,x,concatenation_production_turbines[2,scenario_max_inflow,:],label ="Lower reservoir")
    T

    #Profit
    P=plot(x,concatenation_profit_turbines[1,scenario_max_inflow,:],size=(1200,700),xlabel="N Steps", ylabel="Profit [€]", title="Scenario with maximum inflow",label="Upper reservoir")
    plot!(P,x,concatenation_profit_turbines[2,scenario_max_inflow,:],label="Lower reservoir")
    P
   
    # MINIMUM INFLOW SCENARIO : RESERVOIRS' VOLUMES
    d=plot(x,concatenation_reservoir[1,scenario_min_inflow,:],size=(1200,700),xlabel="N steps",ylabel="Reservoir's volume [Mm3]",title="Scenario with minimum inflow",label = "Upper reservoir")
    plot!(d,x,concatenation_reservoir[2,scenario_min_inflow,:],label = "Lower reservoir")
    d

    # Water Level Variations minimum inflow
    f=plot(x,concatenation_waterLevel_variations[1,scenario_min_inflow,:],size=(1000,500),xlabel="N steps",ylabel="Reservoir's water level variations [Mm3]",title="Scenario with minimum inflow",label="Upper reservoir")
    plot!(f,x,concatenation_waterLevel_variations[2,scenario_min_inflow,:],label="Lower Reservoir")

    # Production
    t=plot(x,concatenation_production_turbines[1,scenario_min_inflow,:],size=(1200,700),xlabel="N Steps", ylabel="Production [MW]", title="Scenario with minimum inflow",label="Upper reservoir")
    plot!(t,x,concatenation_production_turbines[2,scenario_min_inflow,:],label ="Lower reservoir")
    t

    #Profit
    p=plot(x,concatenation_profit_turbines[1,scenario_max_inflow,:],size=(1200,700),xlabel="N Steps", ylabel="Profit [€]", title="Scenario with minimum inflow",label="Upper reservoir")
    plot!(p,x,concatenation_profit_turbines[2,scenario_min_inflow,:],label="Lower reservoir")
    p
    
    
    # SAVING PLOTS
    savefig(C,"Maximum inflow scenario_$scenario_max_inflow.png")
    savefig(c,"Minimum inflow scenario_$scenario_min_inflow.png")
    savefig(D,"Reservoir volumes_$scenario_max_inflow.png")
    savefig(F,"WaterLevel_$scenario_max_inflow.png")
    savefig(T,"Net Production_$scenario_max_inflow.png")
    savefig(P,"Profit scenario_$scenario_max_inflow.png")
    savefig(d,"Reservoir volumes_$scenario_min_inflow.png")
    savefig(f,"WaterLevel_$scenario_min_inflow.png")
    savefig(t,"Net Production_$scenario_min_inflow.png")
    savefig(p,"Profit scenario_$scenario_min_inflow.png")
    
    return("Saved_plots")

end

