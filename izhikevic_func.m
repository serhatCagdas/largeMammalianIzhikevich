function [uout,vout,is_fired] = izhikevic_func(nrn,I,tau,layer_ind, dend_ind)
    
    % the membrane potential of compartment v and recovery variable u are calculated 
    % using comp params 
     
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
    
    dv = ((vout-nrn.vr)*(vout-nrn.vt) - uout + I)*tau/nrn.C ;
    vout = vout + dv;
    
    du = nrn.a*(nrn.b*(vout-nrn.vr) - uout)*tau; 
    uout = uout + du;
     
end

