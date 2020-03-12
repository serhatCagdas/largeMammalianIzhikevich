function [uout,vout,is_fired] = izhikevic_func(nrn,I,tau,layer_ind, dend_ind)
    
    % the membrane potential of compartment v and recovery variable u are calculated 
    % using comp params 
    
    %     exceptional rules
    % 
    %     For LS neurons (nb1)          --> Dendrites are passive 	dv/dt = (-Idend-Isyn) / 250
    %     For LTS neurons (other nb) 	--> always u<=670
    %     For TI 						--> always u<= 530
    %     For TC						--> v>-65 ==> b = 0
    %     For TRN                       --> v>-65 ==> b = 2 
    
    % RS = 1, LS = 2, FS=3, LTS = 4, TS= 5, TI= 6, TRN =7 
    
    is_fired = 0;
    if(dend_ind == 0)
        v_peak = nrn.vpk_soma;
        v = nrn.v;
        u = nrn.u;
        c = nrn.c_soma; 
    else 
        v_peak = nrn.vpk_dend;
        v = nrn.layer{layer_ind}.dend{dend_ind}.v;
        u = nrn.layer{layer_ind}.dend{dend_ind}.u;
        c = nrn.c_dend;
    end
    
    if(v >= v_peak)
        vout = c;
        uout = u + nrn.d; 
        is_fired = 1;
    else
        vout = v;
        uout = u; 
    end    
    
    if (nrn.TypeId ==2 && dend_ind > 0) % exception for nb1 (LS type)
        dv = (I/250)*tau; 
    else                                        % for others 
        dv = ((vout-nrn.vr)*(vout-nrn.vt) - uout + I)*tau/nrn.C ;
    end
    vout = vout + dv;
    
    if (nrn.TypeId == 5 && vout > -65)               %for TC
        du = nrn.a*(-uout)*tau;                      % b = 0
    elseif (nrn.TypeId == 7  && vout > -65)          %for TRN
        du = nrn.a*(2*(vout-nrn.vr) - uout)*tau;     % b = 2
    else
        du = nrn.a*(nrn.b*(vout-nrn.vr) - uout)*tau; % other cases
    end
        uout = uout + du;
    
    if ((nrn.TypeId == 4) && (uout > 670))  % for LTS
        uout = 670;
    elseif ((nrn.TypeId == 6) && (uout > 530)) % for TI
        uout = 530;
    end
     
end

