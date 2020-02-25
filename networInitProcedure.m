function [nrn, smax, isexcitatory ] = networInitProcedure(updateData,fname,tabCellPerc,tabSynapsPerc,tabParam, swParam)

    if(updateData == 0)
        
        load(fname);
        
    else
        disp('initializing network...');
        nrn = cell(0); 
        smax       = []; 
        isexcitatory = [];
        typeId     = []; %need for STD/F
        cellLayer  = [];
        
        for i= 1:length(tabCellPerc.Numb)

            numberOfNeuron = tabCellPerc.Numb(i);
            nrntemp =cell(numberOfNeuron,1);
            smax       = [smax ; ones(numberOfNeuron,1)*tabParam{'smax',tabCellPerc.Row(i)}]; 
            isexcitatory = [isexcitatory ; ones(numberOfNeuron,1)*tabParam{'excitatory',tabCellPerc.Row(i)}];
            typeId       = [typeId ;       ones(numberOfNeuron,1)*tabParam{'typeId',tabCellPerc.Row(i)} ];
            
            for j = 1:numberOfNeuron

                nrntemp{j}.Type = char(tabCellPerc.Row(i));
                nrntemp{j}.C = tabParam{'C',nrntemp{j}.Type}*(1+0.1*randn());   %10 perc rand
                nrntemp{j}.k = tabParam{'k',nrntemp{j}.Type};
                nrntemp{j}.vr = tabParam{'vr',nrntemp{j}.Type}*(1+0.1*rand());
                nrntemp{j}.vt = tabParam{'vt',nrntemp{j}.Type}*(1+0.1*rand());
                nrntemp{j}.vpk_soma = tabParam{'vpk_soma',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.vpk_dend = tabParam{'vpk_dend',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.Gup = tabParam{'Gup',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.Gdown = tabParam{'Gdown',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.a = tabParam{'a',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.b = tabParam{'b',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.c_soma = tabParam{'c_soma',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.c_dend = tabParam{'c_dend',nrntemp{j}.Type}*(1+0.1*randn());
                nrntemp{j}.d = tabParam{'d',nrntemp{j}.Type}*(1+0.1*randn()); 
                
                switch(i)
                    case 1  
                        nrntemp{j}.somaLayer = 1;
                    case  {2,3,4}
                        nrntemp{j}.somaLayer = 23;
                    case  {5,6,7,8,9}
                        nrntemp{j}.somaLayer = 4;
                    case  {10,11,12,13}
                        nrntemp{j}.somaLayer = 5;
                    case  {14,15,16,17}
                        nrntemp{j}.somaLayer = 6;
                    case  {18,20}
                        nrntemp{j}.somaLayer = 8;   %spesific thalamic cells
                    case  {19,21}
                        nrntemp{j}.somaLayer = 9;   %non-spesific thalamic cells
                    case  {22}
                        nrntemp{j}.somaLayer = 10;  %TRN
                end 
                cellLayer = [cellLayer nrntemp{j}.somaLayer];

            end


            nrn = [nrn; nrntemp]; 
        end

        clear numberOfNeuron nrntemp 

        disp('dendrite and synapse establishment...');

        for i= 1:length(nrn)  
                
             if(mod(i,10) == 0)
                 fprintf ('%d\n' ,i);
             else
                 fprintf ('%d ' ,i);
             end
             
             nrn{i}.isfired = 0;
             nrn{i}.v = nrn{i}.c_soma -10;
             nrn{i}.u = nrn{i}.b*nrn{i}.v; 
             nrn{i}.dendLayer  = cell2mat( tabSynapsPerc{nrn{i}.Type  , 'dendLayer'});
             nrn{i}.synapseNum = round(cell2mat (tabSynapsPerc{nrn{i}.Type  ,'synapseNum'}) * swParam.scale_factor *(1 + 0.1*randn()));
            
             nrn{i}.dendNum =  ceil(nrn{i}.synapseNum/swParam.max_synapse_per_dendrite); 
             synapsePerc = tabSynapsPerc(nrn{i}.Type  , :);
             
             for layerInd = 1: length(nrn{i}.synapseNum)
                 
                synapses = zeros(nrn{i}.synapseNum(layerInd),1);
                indPreSyn = 1;
                
                for preSynTypeInd = 1:length(tabCellPerc.Row)
                   
                    firstInd    = tabCellPerc{tabCellPerc.Row(preSynTypeInd),'FirstInd'};
                    cellNum     = tabCellPerc{tabCellPerc.Row(preSynTypeInd),'Numb'}; 
                    synapsePercAll  = cell2mat(synapsePerc{nrn{i}.Type,tabCellPerc.Row(1)});
                    synapseNumLayer = round(synapsePercAll(layerInd)*nrn{i}.synapseNum(layerInd)/100);
                    
                    synapses(indPreSyn : indPreSyn + synapseNumLayer -1 ) ...
                         = firstInd + ceil(rand(synapseNumLayer,1)*cellNum);
                    indPreSyn = indPreSyn + synapseNumLayer;
                end
                
                
                synapseInds = ceil(rand(swParam.max_synapse_per_dendrite,nrn{i}.dendNum(layerInd))*nrn{i}.synapseNum(layerInd));
                
                for dendInd = 1:nrn{i}.dendNum(layerInd)
                    
                    remaining_syn = nrn{i}.synapseNum(layerInd) - swParam.max_synapse_per_dendrite*(dendInd-1);
                     
                    ind_limit = (remaining_syn>swParam.max_synapse_per_dendrite)*swParam.max_synapse_per_dendrite +...
                                (remaining_syn<=swParam.max_synapse_per_dendrite)*remaining_syn;
                    nrn{i}.layer{layerInd}.dend{dendInd}.synapses = ...
                        synapseInds(1:ind_limit ,dendInd);
                     
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.W = ...                           %weights
                        (isexcitatory(nrn{i}.layer{layerInd}.dend{dendInd}.synapses)).*rand(ind_limit,1)*6 ... % for excitatory
                        +(1 - isexcitatory(nrn{i}.layer{layerInd}.dend{dendInd}.synapses))*4; % for inhibitory
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.gNMDA = nrn{i}.layer{layerInd}.dend{dendInd}.W .* ... 
                                                                 (isexcitatory(nrn{i}.layer{layerInd}.dend{dendInd}.synapses));
                    nrn{i}.layer{layerInd}.dend{dendInd}.gAMPA = nrn{i}.layer{layerInd}.dend{dendInd}.gNMDA;
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.gGABAa = nrn{i}.layer{layerInd}.dend{dendInd}.W .* ...  
                                                                 (1 - isexcitatory(nrn{i}.layer{layerInd}.dend{dendInd}.synapses));
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.gGABAb = nrn{i}.layer{layerInd}.dend{dendInd}.gGABAa ;
                    
                    [nrn{i}.layer{layerInd}.dend{dendInd}.STD.Ts,...
                        nrn{i}.layer{layerInd}.dend{dendInd}.STD.p] ...
                        = genShortTimeParamVector(typeId(i),typeId(nrn{i}.layer{layerInd}.dend{dendInd}.synapses));
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.STD.x = ones(ind_limit,1);
                    nrn{i}.layer{layerInd}.dend{dendInd}.STDP.c     = zeros(ind_limit,1);
                    nrn{i}.layer{layerInd}.dend{dendInd}.fired_ind  = -1000;
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.v = nrn{i}.c_dend-10;
                    nrn{i}.layer{layerInd}.dend{dendInd}.u = nrn{i}.c_dend*nrn{i}.b;
                    
                    nrn{i}.layer{layerInd}.dend{dendInd}.delay = calcAxonalDelay(nrn{i}.dendLayer(layerInd) ,cellLayer(nrn{i}.layer{layerInd}.dend{dendInd}.synapses));
                    
                end 
                    
             end
            
        end
        
        description = {'10 perc randomization added','different inital v'};
        
        save(fname,'nrn','smax','isexcitatory','description');
    end

end

