function [tabCellPerc,tabSynapsPerc,tabParam] = generateTables(swParam)


    cellType = {'nb1','p23','b23','nb23','ss4L4','ss4L23','p4','b4','nb4','p5L23','p5L56','b5','nb5','p6L4',...
                'p6L56','b6','nb6','TCs','TCn','TIs','TIn','TRN'}; 
    cellPerc = [1.5, 26, 3.1, 4.2, 9.2, 9.2, 9.2, 5.4, 1.5, 4.8, 1.3, 0.6, 0.8, 13.6,...
                4.5, 2, 2, 0.5, 0.5, 0.1, 0.1, 0.5];        
    cellNumb = floor(cellPerc* swParam.Nnrn/100);
    cumulCellNumb = cumsum(cellNumb);
    Indices = [1 cumulCellNumb(1:end-1)+1];
    tabCellPerc = table(cellPerc',cellNumb',Indices');
    tabCellPerc.Properties.RowNames = cellType; 
    tabCellPerc.Properties.VariableNames{'Var1'} = 'Perc';
    tabCellPerc.Properties.VariableNames{'Var2'} = 'Numb';
    tabCellPerc.Properties.VariableNames{'Var3'} = 'FirstInd';

    load('synapseArray.mat');
    dendLayer = {1,[23,1],23,23,4,4,[4,23,1],4,4,[5,4,23,1],[5,4,23,1],5,5,[6,5,4,23],[6,5,4,23,1],6,6,8,9,8,9,10}; 
    newArr = [dendLayer' synapseArray];

    tabSynapsPerc = cell2table(newArr,...
        'VariableNames',['dendLayer', 'synapseNum' ,cellType(1:17),'corCort',cellType(18:22),'brainstem','sensory']);
    tabSynapsPerc.Properties.RowNames = cellType;

    load('paramArray.mat');
    paramName ={'C','k','vr','vt','vpk_soma','vpk_dend','Gup','Gdown','a','b','c_soma','c_dend','d','excitatory','smax','typeId'};
    % newArr = [paramName' paramArray];  
    isExcitatory = {0 1  0 0  1  1  1 0 0  1  1 0 0  1  1 0 0  1  1 0 0 0};
    smax         = {4 10 6 4 10 10 10 6 4 10 10 6 4 10 10 6 4 20 20 5 5 5}; 
    %p:1,ss:2,b:3,nb:4,TC:5,TI:6,RTN:7
    typeId = {4,1,3,4,2,2,1,3,4,1,1,3,4,1,1,3,4,5,5,6,6,7};
    
    tabParam = cell2table([paramArray ; isExcitatory; smax; typeId],...
        'VariableNames', cellType);
    tabParam.Properties.RowNames = paramName; 
    
end

