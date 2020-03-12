function [] = print_raster(data,title_d, T,tau, isplot,issave)
 
    MarkerFormat.MarkerSize = 4;
    MarkerFormat.Color = [0 0 1];
    figure, 
    plotSpikeRaster(data ,'PlotType','scatter','XLimForCell',[0 T],'TimePerBin',tau,'MarkerFormat',MarkerFormat);
    ylabel('Trial') 
    title(title_d, 'Interpreter', 'none')


    if(issave)
        saveas(gcf,title_d,'png')
    end
    if(isplot~=1)
        close all;
    end
         
    
end

