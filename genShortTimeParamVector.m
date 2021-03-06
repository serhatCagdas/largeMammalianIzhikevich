function [Ts,p] = genShortTimeParamVector( postSynNrn ,preSynNrns)

synType = 10*preSynNrns + postSynNrn ;

% RS = 1, LS = 2, FS=3, LTS = 4, TS= 5, TI= 6, TRN =7 

 state150_06  = (synType == 11)  + (synType == 13) ...
              + (synType == 31) + (synType == 33);

 state150_07 = (synType == 51) ;

 state100_15 = (synType == 12) + (synType == 14);
 
 state200_06 = (synType == 53);
 
 rest = 1 - (state150_06 + state150_07 + state100_15 + state200_06);
 
 Ts = (state150_06+state150_07)*150 + state100_15*100 + state200_06*200 + rest;
 
 Ts = Ts.*(1 + 0.1*rand(length(Ts),1));
 
 p  = (state150_06 + state200_06)*0.6 + state150_07*0.7 + state100_15*1.5 + rest;
 p  = p.*(1 + 0.1*randn(length(p),1));
 
end

