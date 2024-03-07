function plot_colorbar(options)
% this is the function to generate colorbar without tick labels
% the input will be only options
    arguments
        options.cmap = 'jet'
        options.size (1,2) double = [1000 100]
        options.vertical = false
    end
    % hf = figure('Units','normalized'); 
    % colormap(options.cmap)
    % if ~options.vertical
    %     hCB = colorbar('north','ticks',[]);
    %     set(gca,'Visible',false)
    %     hCB.Position = options.position;
    %     hf.Position(4) = 0.1000;
    % else
    %     hCB = colorbar('east','ticks',[]);
    %     set(gca,'Visible',false)
    %     hCB.Position = [0.3 0.1 0.4 0.8];
    %     hf.Position(3) = 0.1000;
    % end

    hf = figure();
    hf.Position(3:4) = options.size + [100 100];
    if options.size(1) > options.size(2)
        hCB = colorbar('north','Ticks',[]);
    else
        hCB = colorbar('east','Ticks',[]);
    end
    colormap(options.cmap);
    
    set(gca,'Visible',false)
    set(hCB,"Units",'pixels',"Position",[50 50 options.size(1) options.size(2)], ...
        "LineWidth",2)

end