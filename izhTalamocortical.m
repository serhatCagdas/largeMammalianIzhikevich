clc;clear all;close all;
%% TASARIM ASAMALARI
% Ref. : Izhikevich, Eugene M., and Gerald M. Edelman. "Large-scale model 
% of mammalian thalamocortical systems." Proceedings of the national academy 
% of sciences 105.9 (2008): 3593-3598.
 
%% parameters 

T = 100; %(ms): 1 sec simulation
Tau.AMPA  = 5;  Tau.NMDA  = 150;
Tau.GABAa = 6;  Tau.GABAb = 150;
t_step = 0.125; % 1/8 ms
A_stdp_p    = 1;    A_stdp_n    = -2;
Tau_stdp_p  = 20;   Tau_stdp_n  = 20; %(ms)
Tau_dop = 200;      Tau_stdp_c = 1000; %(ms)
DA_tonic = 0.01;
d = 0.002;      %concentration of extracellular DA
                % baseline concentration

N = round(T*1/t_step);

Nreward = ceil(T/500); % average of 1 reward per 500 ms
n_reward = ceil(rand(Nreward,1)*N); 

swParam.Nnrn = 1000; 
swParam.max_synapse_per_dendrite =  40;
swParam.scale_factor             =  0.05;
updateInit = 1;

[tabCellPerc,tabSynapsPerc,tabParam] = generateTables(swParam);

fname = 'nrn1k_200224_del1.mat';
[nrn , smax, isexcitatory] = networInitProcedure(updateInit, fname ,tabCellPerc,tabSynapsPerc,tabParam, swParam);
clear fname updateInit

disp('simulation started...');

Nnrn = length(nrn);
spikingTimes=cell(Nnrn,1);
firedtmp = zeros(Nnrn,1);
fired      = firedtmp;
fired_inds = -1000*ones(Nnrn,1);
r.AMPA  = ones(Nnrn,1); 
r.NMDA  = ones(Nnrn,1);
r.GABAa = ones(Nnrn,1);
r.GABAb = ones(Nnrn,1);

AL = 20/t_step; %axon length 20 ms max
axon = logical(zeros(Nnrn,AL));

x = ones(Nnrn,1); % STD/F parameter

% STDP = A_stdp_n*(t<0).*exp(t/Tau_stdp_n) + A_stdp_p*(t>=0).*exp(-t/Tau_stdp_p);

for n = 1:N 
    
    %Dopamine regulation
    DA =  0.01 + sum(n == n_reward)*((0.5*t_step)-0.01); % if reward DA = 0.5 else 0.01 
    d = d + t_step*(-d/Tau_dop + DA);
    
    % network dynamic
    for nind = 1:Nnrn
%          nind
        nrntmp = nrn{nind};
        
        %calc v u for soma 
        Idend = 0;
        Isyn  = 0;
        for d_ind = 1:length(nrn{nind}.layer{1})
            Idend = Idend + nrntmp.Gup*(nrntmp.v - nrn{nind}.layer{1}.dend{d_ind}.v );
        end 
        I = -Idend-Isyn;
        [nrntmp.u,nrntmp.v,is_fired] = izhikevic_func(nrn{nind},I,t_step,1, 0);
        
        if(is_fired) spikingTimes{nind} = [spikingTimes{nind} n*t_step]; end
        firedtmp(nind) = is_fired;
        
        % v u for dendritic comps
        for l_ind = 1:length(nrn{nind}.dendLayer) %each layer
%              l_ind
            for d_ind = 1:length(nrn{nind}.layer{l_ind})
%                  d_ind
                Idend = 0;
                Isyn = 0; 
                 
                % calc dendritic currents caused of mother compartments 
                if(l_ind == 1) 
                    Idend = Idend + nrntmp.Gdown*( nrn{nind}.layer{1}.dend{d_ind}.v - nrn{nind}.v);
                else
                    Idend = Idend + nrntmp.Gdown*( nrn{nind}.layer{1}.dend{d_ind}.v - ... %former layer first dend is mother
                                            nrn{nind}.layer{l_ind-1}.dend{1}.v);
                end
                
                % calc dendritic currents caused of daugther compartments 
                if (d_ind == 1 && l_ind < length(nrn{nind}.dendLayer)) 
                    
                    for d_ind_up = 1:length(nrn{nind}.layer{l_ind+1})
                        Idend = Idend + nrntmp.Gup*(nrn{nind}.layer{l_ind}.dend{1}.v - ...
                                            nrn{nind}.layer{l_ind+1}.dend{d_ind_up}.v );
                    end 
                end
                
                % calc synaptic current
                dend = nrn{nind}.layer{l_ind}.dend{d_ind};
                
                fired = (axon(dend.synapses' + (Nnrn * floor(nrn{nind}.layer{l_ind}.dend{d_ind}.delay/t_step))))';
                 
                Isyn =sum( (fired .* dend.STD.x).*...
                          ( r.AMPA(dend.synapses).*dend.gAMPA*(dend.v) + ...
                            r.NMDA(dend.synapses).*dend.gNMDA*(((dend.v+80)/60)^2)/...
                                                               (1+((dend.v+80)/60)^2)* dend.v+...
                            r.GABAa(dend.synapses).*dend.gGABAa*(dend.v+70) +...
                            r.GABAb(dend.synapses).*dend.gGABAb*(dend.v+90)));
                                               
                
                I = -Idend-Isyn;
                [uout,vout,is_fired] = izhikevic_func(nrn{nind},I,t_step,l_ind, d_ind);
                nrntmp.layer{l_ind}.dend{d_ind}.v = vout;
                nrntmp.layer{l_ind}.dend{d_ind}.u = vout; 
                
                % STD/F : Short time depression / facilitation
                dend.STD.x = fired.* dend.STD.x .* dend.STD.p + ...      %if fired
                            (1- fired) .* (dend.STD.x + t_step*(1-dend.STD.x)./dend.STD.Ts);  %if not
                
                        
                %STDP : Spike timing dependent plasticity
%                 dend.STDP.count = (1- fired(dend.synapses)).*(dend.STDP.count + t_step); % if fired set to zero
                
                
                if(is_fired)
                    dend.fired_ind = n;
%                     STDP = A_stdp_n*(dend.STDP.count<0).*exp(dend.STDP.count/Tau_stdp_n) +...
%                            A_stdp_p*(dend.STDP.count>=0).*exp(-dend.STDP.count/Tau_stdp_p);
                    diff_t = t_step*(fired_inds(dend.synapses)-n);
                    STDP = A_stdp_p *exp(diff_t/Tau_stdp_p);
                else
                    
                    diff_t = t_step*(dend.fired_ind-fired_inds(dend.synapses));
                    STDP = -fired.*(A_stdp_n * exp(diff_t/Tau_stdp_n));
                end
                 
                dend.STDP.c = dend.STDP.c * (1- t_step/Tau_stdp_c)+ STDP; 
                dend.STDP.sd = d * dend.STDP.c; 
                dend_smax = smax(dend.synapses);
                
                %update conductance
                dend.gAMPA = dend.gAMPA + (dend.gAMPA < dend_smax).*(t_step*dend.STDP.sd); 
                dend.gNMDA = dend.gNMDA + (dend.gNMDA < dend_smax).*(t_step*dend.STDP.sd);  
                
            end
        end
        
        nrn{nind} = nrntmp;
    end
    if(mod(n,10) == 0)
         fprintf ('%d\n' ,n);
     else
         fprintf ('%d ' ,n);
     end
    
    fired = firedtmp;
    axon = [fired axon(:,1:AL-1)];
    fired_inds = (fired==1)*n + (fired==0).*fired_inds;
    
    r.AMPA = (fired.*ones(Nnrn,1)) + ((1-fired).* r.AMPA *(1 - t_step/Tau.AMPA));
    r.NMDA = (fired.*ones(Nnrn,1)) + ((1-fired).* r.NMDA *(1 - t_step/Tau.NMDA));
    r.GABAa = (fired.*ones(Nnrn,1)) + ((1-fired).*r.GABAa*(1 - t_step/Tau.GABAa));
    r.GABAb = (fired.*ones(Nnrn,1)) + ((1-fired).*r.GABAb*(1 - t_step/Tau.GABAb));
    
%     x = 
     
end


for nind = 1:Nnrn
   
    if(isempty(spikingTimes{nind}))
        spikingTimes{nind} = 0;
    end 
end

print_raster(spikingTimes,'', T,t_step, 1,0);













