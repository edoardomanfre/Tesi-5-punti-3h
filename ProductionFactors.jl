# EVALUATION OF PRODUCTION FACTORS AND PRODUCTION TIME

using Plots
function production_factors(InputParameters::InputParam)

    max_step_production_turbines=zeros(HY.NMod);   
    max_step_request_pump=0;                                            
    NStep=InputParameters.NStep

    for i=1:HY.NMod
        a=findmax(HY.Eff[i,:])                                                                                                 # Per ogni turbina, trovo la massima efficienza e la corrispettiva Potenza MW prodotta
        max_step_production_turbines[i]=HY.Eff[i,a[2]]*HY.DisMaxSeg[i,a[2]]
    end

    b=findmax(HY.EffPump[:])
    max_step_request_pump=1/HY.EffPump[b[2]]*HY.DisMaxSegPump[b[2]]
    

    step_production_factor_turbines=zeros(HY.NMod,NSimScen,NStage,NStep);
    step_factor_pump=zeros(HY.NMod,NSimScen,NStage,NStep);
    
    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            for iStage=1:NStage
                for iStep=1:NStep
                step_production_factor_turbines[iMod,iScen,iStage,iStep]=ResultsSim.Production[iMod,iScen,iStage,iStep]/max_step_production_turbines[iMod]                         # fattore di produzione giornaliero per ogni turbina, ogni scenario
                step_factor_pump[iMod,iScen,iStage,iStep]=ResultsSim.Pumping[iMod,iScen,iStage,iStep]/max_step_request_pump
                end
                #weekly_production_factor[i,s,t]=weekly_production[i,s,t]/max_weekly_production[i]                               #Fattore di produzione settimanale
            end
        end
    end


    # CONCATENATION OF PRODUCTION AND PRODUCTION FACTOR
    concatenation_production_turbines=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_production_pump=zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_productionFactor_turbines= zeros(HY.NMod,NSimScen,NStep*NStage);
    concatenation_productionFactor_pumps=zeros(HY.NMod,NSimScen,NStep*NStage);

    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            for iStage=1:NStage
                concatenation_production_turbines[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.Production[iMod,iScen,iStage,:]
                concatenation_production_pump[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=ResultsSim.Pumping[iMod,iScen,iStage,:]
                concatenation_productionFactor_turbines[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=step_production_factor_turbines[iMod,iScen,iStage,:]
                concatenation_productionFactor_pumps[iMod,iScen,(NStep*(iStage-1))+1:1:(NStep*iStage)]=step_factor_pump[iMod,iScen,iStage,:]
            end
        end
    end

    #How many days, in a given scenario, have a certain capacity factor (Production_step=0.25)

    nsteps_turbines=zeros(HY.NMod,NSimScen,7);
    nsteps_pumps=zeros(HY.NMod,NSimScen,7);
    count_scenarios=zeros(HY.NMod,2);
    for iMod=1:HY.NMod
        for iScen=1:NSimScen
                nsteps_turbines[iMod,iScen,1]=count(i->(i==0),step_production_factor_turbines[iMod,iScen,:,:])  
                nsteps_turbines[iMod,iScen,2]=count(i->(0<i<0.25),step_production_factor_turbines[iMod,iScen,:,:])                                                             #How many steps work for P< 0.25
                nsteps_turbines[iMod,iScen,3]=count(i->(i<=0.5 && i>=0.25),step_production_factor_turbines[iMod,iScen,:,:])
                nsteps_turbines[iMod,iScen,4]=count(i->(i<=0.75 && i>0.5),step_production_factor_turbines[iMod,iScen,:,:])
                nsteps_turbines[iMod,iScen,5]=count(i->(i<1.00 && i >0.75),step_production_factor_turbines[iMod,iScen,:,:])
                nsteps_turbines[iMod,iScen,6]=count(i->(i==1.00),step_production_factor_turbines[iMod,iScen,:,:])
                nsteps_turbines[iMod,iScen,7]=count(i->(i>1.00),step_production_factor_turbines[iMod,iScen,:,:])

                nsteps_pumps[iMod,iScen,1]=count(i->(i==0),step_factor_pump[iMod,iScen,:,:])  
                nsteps_pumps[iMod,iScen,2]=count(i->(0<i<0.25),step_factor_pump[iMod,iScen,:,:])                                                             #How many steps (7 or 56) work for P< 0.25
                nsteps_pumps[iMod,iScen,3]=count(i->(i<=0.5 && i>=0.25),step_factor_pump[iMod,iScen,:,:])
                nsteps_pumps[iMod,iScen,4]=count(i->(i<=0.75 && i>0.5),step_factor_pump[iMod,iScen,:,:])
                nsteps_pumps[iMod,iScen,5]=count(i->(i<1.00 && i >0.75),step_factor_pump[iMod,iScen,:,:])
                nsteps_pumps[iMod,iScen,6]=count(i->(i==1.00),step_factor_pump[iMod,iScen,:,:])
                nsteps_pumps[iMod,iScen,7]=count(i->(i>1.00),step_factor_pump[iMod,iScen,:,:])
        end
    end
    nsteps_turbines;
    nsteps_pumps;

    # Calculates the n. of scenarios which are working at MAximum Efficiency Point and ones who are not working at all for more than 50% of time in the year

    for iMod=1:HY.NMod
        count_scenarios[iMod,1]= count(i->(i>=0.5),nsteps_turbines[iMod,:,6]/(NStep*NStage))                                      # Ones working at Max Efficiency
        count_scenarios[iMod,2]= count(i->(i>=0.5),nsteps_turbines[iMod,:,1]/(NStep*NStage))                                      # Those which are not working at all
    end

    #Do the mean for all scenarios
    mean_nsteps_turbines=zeros(HY.NMod,7);
    mean_nsteps_pumps=zeros(HY.NMod,7);
    for iMod=1:HY.NMod
        for i=1:7
            mean_nsteps_turbines[iMod,i]=mean(nsteps_turbines[iMod,:,i])
            mean_nsteps_pumps[iMod,i]=mean(nsteps_pumps[iMod,:,i])
        end
    end

    b=zeros(HY.NMod,NSimScen,NStage*NStep);
    c=zeros(HY.NMod,NSimScen,NStage*NStep);  
    for iMod=1:HY.NMod
        for iScen=1:NSimScen
            b[iMod,iScen,:]=sort(concatenation_productionFactor_turbines[iMod,iScen,:],rev=true)
            c[iMod,iScen,:]=sort(concatenation_productionFactor_pumps[iMod,iScen,:],rev=true)
        end
    end

    return ProductionTime(
        max_step_production_turbines,
        max_step_request_pump,
        step_production_factor_turbines,
        step_factor_pump,
        nsteps_turbines,
        nsteps_pumps,
        mean_nsteps_turbines,
        mean_nsteps_pumps,
       # best_scenarios_turbines,
        #worst_scenarios_turbines,
        #best_scenarios_pumps,
        #worst_scenarios_pumps,
        #count_scenarios,
        #profit,
        #plot_turbines,
        #plot_pumps,
        )
end

  