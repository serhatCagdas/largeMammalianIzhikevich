%  in gray  matter (non-mylinated)                0.1 m/s = 0.1 mm/ms
%  in white matter (mylinated)                    1   m/s
%  cortex is 1 mm witdth  from layer 6 to 1   =>  10  ms (since non myelinated)
%  + L5 to thalamus                               1   ms
%  + L6 to thalamus                               20  ms
%  + specific thalamocortical                     1   ms
%  - non-specific thalamocortical determined according to axonal length (no
%  longer than 20 ms) (randomized for now)


function [delay_ms] = calcAxonalDelay(postDendLayer,preSomaLayer)

    postDendLayer(postDendLayer == 23) = 2.5;
    preSomaLayer (preSomaLayer  == 23) = 2.5;
    
    if( postDendLayer == 8 || postDendLayer == 9)
        
        delay_ms = (preSomaLayer == 5) + (preSomaLayer == 6)*20; %  + L5 to thalamus	1  ms
                                                                 %  + L6 to thalamus	20 ms
        
        delay_ms = delay_ms + (postDendLayer == 8); %  specific thalamocortical 	1 ms
        
        delay_ms = delay_ms + (postDendLayer == 9).*rand(1,length(delay_ms))*20; %  non specific randomized for now
        
    else
        
        delay_ms = abs(postDendLayer - preSomaLayer).*(1+0.1*randn(1,length(preSomaLayer)));        
    end

end

