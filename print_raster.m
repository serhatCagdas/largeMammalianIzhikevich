function [] = print_raster(data,title_d, T,tau, isplot,issave)
 
    figure, 
    plotSpikeRaster(data ,'PlotType','vertline','XLimForCell',[0 T],'TimePerBin',tau);
    ylabel('Trial') 
    title(title_d, 'Interpreter', 'none')


    if(issave)
        saveas(gcf,title_d,'png')
    end
    if(isplot~=1)
        close all;
    end
         
    
end

