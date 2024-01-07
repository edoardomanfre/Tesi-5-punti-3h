# DEPENDENCE ON DOWNSTREAM WATER LEVEL

function DeactivationPump_sim(SP,iScen,t,HY,Reservoir,limit,NStep)

    reservoir = 0

    if t==1 
        if iScen==1
            if HY.ResInit0[2] <= limit  #52
                for iStep=1:NStep
                    #JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],0)     
                    JuMP.set_normalized_rhs(SP.maxReleasePump[iStep],0)
                    reservoir = HY.ResInit0[2]
                    #SP = add_disLimitPump(SP,NStep)
                end
            elseif HY.ResInit0[2] > limit
                for iStep=1:NStep
                    #JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],HY.DisMaxSegPump[tSeg]) 
#                    JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], sum(HY.DisMaxSegPump[tSeg] for tSeg = 1:HY.NDSegPump))
                    JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], HY.DisPointPump[2])
                    #SP= relax_disLimitPump(SP,NStep)
                    reservoir = HY.ResInit0[2]
                end
            end
        else    #iScen>1
            if Reservoir[2,iScen-1,end,end] <= limit
                for iStep=1:NStep
                    #JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],0)     
                    JuMP.set_normalized_rhs(SP.maxReleasePump[iStep],0)
                    reservoir = Reservoir[2,iScen-1,end,end]
                    #SP = add_disLimitPump(SP,NStep)
                end
            elseif Reservoir[2,iScen-1,end,end] > limit
                for iStep=1:NStep
                    #JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],HY.DisMaxSegPump[tSeg]) 
#                    JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], sum(HY.DisMaxSegPump[tSeg] for tSeg = 1:HY.NDSegPump))
                    JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], HY.DisPointPump[2])
                    reservoir = Reservoir[2,iScen-1,end,end]
                    #SP= relax_disLimitPump(SP,NStep)
                end
            end
        end
    else    #if t>1 50 scen4
        if Reservoir[2,iScen,t-1,end] <= limit
           for iStep=1:NStep
                #JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],0)
                JuMP.set_normalized_rhs(SP.maxReleasePump[iStep],0)
                reservoir = Reservoir[2,iScen,t-1,end]
                #SP = add_disLimitPump(SP,NStep)
            end
        elseif Reservoir[2,iScen,t-1,end] > limit
            for iStep=1:NStep
               # JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],HY.DisMaxSegPump[tSeg])
#                JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], sum(HY.DisMaxSegPump[tSeg] for tSeg = 1:HY.NDSegPump))
                JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], HY.DisPointPump[2])
                reservoir = Reservoir[2,iScen,t-1,end]
                #SP= relax_disLimitPump(SP,NStep)
            end
        end
    end

    return SP,reservoir

end

function DeactivationPump_SDP(SP,HY,ResSeg,limit,nfrom,NStep)         # iMod=2 

        if ResSeg[nfrom][2] <= limit
            for iStep=1:NStep
             #   JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],0)
                JuMP.set_normalized_rhs(SP.maxReleasePump[iStep],0)
                #SP = add_disLimitPump(SP,NStep)
            end
        elseif ResSeg[nfrom][2] > limit
            for iStep=1:NStep
               #JuMP.set_normalized_rhs(SP.pumpdischarge[tSeg,iStep],HY.DisMaxSegPump[tSeg])
                JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], sum(HY.DisMaxSegPump[tSeg] for tSeg = 1:HY.NDSegPump))
#                JuMP.set_normalized_rhs(SP.maxReleasePump[iStep], HY.DisPointPump[2])
                #SP= relax_disLimitPump(SP,NStep)
            end
        end

    return SP

end

function add_disLimitPump(SP, NStep)                                    # attivo i vincoli sullo scarico                                                                         
    for iStep = 1:NStep
      set_normalized_rhs(SP.maxReleasePump[iStep], 0)
    end
    return SP
end

function relax_disLimitPump(SP,NStep)                                   # disattivo i vincoli sullo scarico
#=    for iStep = 1:NStep
      set_normalized_rhs(
        SP.maxReleasePump[iStep],
        sum(HY.DisMaxSegPump[tSeg] for tSeg = 1:HY.NDSegPump),
      )
    end=#
    for iStep = 1:NStep
        set_normalized_rhs(
          SP.maxReleasePump[iStep],
          HY.DisPointPump[2],
        )
      end
      
    return SP
end



# RAMPING CONSTRAINTS - 

function intra_volume_changes(case::caseData,SP,iMod,iScen,t,iStep,HY,Reservoir,NStep)
    
    path=case.DataPath
    cd(path)
    f=open("Water_volumes_levels.dat")
    line=readline(f)

    line = readline(f)
    items = split(line, " ")
    NMod = parse(Int, items[1]) #set number of modules
    NVolumes=zeros(NMod);
 
    Volume_changes= zeros(HY.NMod,Int(NVolumes[iMod])-1)
    Volume_changes[1,:]=[0.0018 0.0042 0.0038 0.0042 0.005 0.005 0.008 0.01 0.01 0.01 0.012 0.014 0.016 0.014 0.012 0.022 0.022 0.022 0.022 0.026]                   
    Volume_changes[2,:]=[0.005 0.00625 0.0055 0.007 0.00875 0.005 0.0065 0.0085 0.0075 0.01 0.01 0.0125]                    

    water_volumes_file=zeros(HY.NMod,Int(NVolumes[iMod]))     
    water_volumes_file[1,:]= [0.00 0.09 0.30 0.49 0.70 0.95 1.20 1.60 2.10 2.60 3.10 3.70 4.40 5.20 5.90 6.50 7.60 8.70 9.80 10.90 12.20]        
    water_volumes_file[2,:]= [0.00 0.20 0.45 0.67 0.95 1.30 1.50 1.76 2.10 2.40 2.80 3.20 3.70]           


    for n=1:size(Volume_changes)[2]
        if t==1
            if HY.ResInit0[iMod]>= water_volumes_file[iMod,n] && HY.ResInit0[iMod]< water_volumes_file[iMod,n+1]
                JuMP.set_normalized_rhs(SP.positive_var[iMod,iStep],Volume_changes[iMod,n])
                JuMP.set_normalized_rhs(SP.negative_var[iMod,iStep],-Volume_changes[iMod,n])
            end
        else #t>1
            if Reservoir[iMod,iScen,t-1,NStep]>= water_volumes_file[iMod,n] && Reservoir[iMod,iScen,t-1,NStep] < water_volumes_file[iMod,n+1]
                JuMP.set_normalized_rhs(SP.positive_var[iMod,iStep],Volume_changes[iMod,n])
                JuMP.set_normalized_rhs(SP.negative_var[iMod,iStep],-Volume_changes[iMod,n])
            end
        end
   end

    return SP
end




function initial_volume_changes(case::caseData,SP,iMod,iScen,t,HY,Reservoir)
    
    path=case.DataPath
    cd(path)
    f=open("Water_volumes_levels.dat")
    line=readline(f)

    line = readline(f)
    items = split(line, " ")
    NMod = parse(Int, items[1]) #set number of modules
    NVolumes=zeros(NMod);
  
    Volume_changes= zeros(HY.NMod,Int(NVolumes[iMod])-1)
    Volume_changes[1,:]=[0.0018 0.0042 0.0038 0.0042 0.005 0.005 0.008 0.01 0.01 0.01 0.012 0.014 0.016 0.014 0.012 0.022 0.022 0.022 0.022 0.026]
    Volume_changes[2,:]=[0.005 0.00625 0.0055 0.007 0.00875 0.005 0.0065 0.0085 0.0075 0.01 0.01 0.0125]  

    water_volumes_file=zeros(HY.NMod,Int(NVolumes[iMod]))
    water_volumes_file[1,:]= [0.00 0.09 0.30 0.49 0.70 0.95 1.20 1.60 2.10 2.60 3.10 3.70 4.40 5.20 5.90 6.50 7.60 8.70 9.80 10.90 12.20] 
    water_volumes_file[2,:]= [0.00 0.20 0.45 0.67 0.95 1.30 1.50 1.76 2.10 2.40 2.80 3.20 3.70]  


    for n=1:size(Volume_changes)[2]       
        if iScen==1
            if t==1
                if HY.ResInit0[iMod]>= water_volumes_file[iMod,n] && HY.ResInit0[iMod] < water_volumes_file[iMod,n+1]
                JuMP.set_normalized_rhs(SP.Initial_volumevar_positive[iMod],HY.ResInit0[iMod]+Volume_changes[iMod,n])
                JuMP.set_normalized_rhs(SP.Initial_volumevar_negative[iMod],HY.ResInit0[iMod]-Volume_changes[iMod,n])
                end
            else # t>1
                if Reservoir[iMod,iScen,t-1,end]>= water_volumes_file[iMod,n] && Reservoir[iMod,iScen,t-1,end] < water_volumes_file[iMod,n+1]
                JuMP.set_normalized_rhs(SP.Initial_volumevar_positive[iMod],Reservoir[iMod,iScen,t-1,end]+Volume_changes[iMod,n] )       
                JuMP.set_normalized_rhs(SP.Initial_volumevar_negative[iMod],Reservoir[iMod,iScen,t-1,end]-Volume_changes[iMod,n] )       
                end
            end
        else            #iScen >1
            if t==1
                if Reservoir[iMod,iScen-1,end,end]>= water_volumes_file[iMod,n] && Reservoir[iMod,iScen-1,end,end] < water_volumes_file[iMod,n+1]
                    JuMP.set_normalized_rhs(SP.Initial_volumevar_positive[iMod],Reservoir[iMod,iScen-1,end,end]+Volume_changes[iMod,n] )       
                    JuMP.set_normalized_rhs(SP.Initial_volumevar_negative[iMod],Reservoir[iMod,iScen-1,end,end]-Volume_changes[iMod,n] )        
                end
            else        #t>1
                if Reservoir[iMod,iScen,t-1,end]>= water_volumes_file[iMod,n] && Reservoir[iMod,iScen,t-1,end] < water_volumes_file[iMod,n+1]
                    JuMP.set_normalized_rhs(SP.Initial_volumevar_positive[iMod],Reservoir[iMod,iScen,t-1,end]+Volume_changes[iMod,n] )          
                    JuMP.set_normalized_rhs(SP.Initial_volumevar_negative[iMod],Reservoir[iMod,iScen,t-1,end]-Volume_changes[iMod,n] )          
                end
            end
        end
    end

    return SP
end


function ramping_constraints_SDP(case::caseData,SP,iMod,t,nfrom,HY,ResSeg,iStep)
    
    path=case.DataPath
    cd(path)
    f=open("Water_volumes_levels.dat")
    line=readline(f)

    line = readline(f)
    items = split(line, " ")
    NMod = parse(Int, items[1]) #set number of modules
    NVolumes=zeros(NMod);
 
    Volume_changes= zeros(HY.NMod,Int(NVolumes[iMod])-1) #3cm variations
    Volume_changes[1,:]=[0.0018 0.0042 0.0038 0.0042 0.005 0.005 0.008 0.01 0.01 0.01 0.012 0.014 0.016 0.014 0.012 0.022 0.022 0.022 0.022 0.026]                           
    Volume_changes[2,:]=[0.005 0.00625 0.0055 0.007 0.00875 0.005 0.0065 0.0085 0.0075 0.01 0.01 0.0125]                 

    water_volumes_file=zeros(HY.NMod,Int(NVolumes[iMod]))
    water_volumes_file[1,:]= [0.00 0.09 0.30 0.49 0.70 0.95 1.20 1.60 2.10 2.60 3.10 3.70 4.40 5.20 5.90 6.50 7.60 8.70 9.80 10.90 12.20]               
    water_volumes_file[2,:]= [0.00 0.20 0.45 0.67 0.95 1.30 1.50 1.76 2.10 2.40 2.80 3.20 3.70]                   


    for n=1:size(Volume_changes)[2]             # iterating from 1 to 4
       
        if iStep==1
            if ResSeg[nfrom][iMod]>= water_volumes_file[iMod,n] && ResSeg[nfrom][iMod] < water_volumes_file[iMod,n+1]
                JuMP.set_normalized_rhs(SP.Initial_volumevar_positive[iMod],ResSeg[nfrom][iMod]+Volume_changes[iMod,n] )       
                JuMP.set_normalized_rhs(SP.Initial_volumevar_negative[iMod],ResSeg[nfrom][iMod]-Volume_changes[iMod,n] )       
            end   
        elseif iStep>1 
            if ResSeg[nfrom][iMod]>= water_volumes_file[iMod,n] && ResSeg[nfrom][iMod] < water_volumes_file[iMod,n+1]
                JuMP.set_normalized_rhs(SP.positive_var[iMod,iStep],Volume_changes[iMod,n])
                JuMP.set_normalized_rhs(SP.negative_var[iMod,iStep],-Volume_changes[iMod,n])
            end
        end

    end
   
    return SP
end