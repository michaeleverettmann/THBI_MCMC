function [] = plot_corrplot(par, allmodels_collated, figname)

    m = allmodels_collated; 
    figure(4082); clf; hold on; hF = gcf; hF.Position(3:4) = [800, 700];
    [~, ~, h] = corrplot([[m.cxi]', [m.mxi]'], 'VarNames', ["xi crust", "xi"+string(num2str(par.mod.mantle.xidepths))']); 
    lineHandles = h(strcmp(get(h, 'type'), 'line'));
    

    %%% The axis labels are weird, by default. Use this code to fix it, https://www.mathworks.com/matlabcentral/answers/172723-corrplot-plotting-on-strange-x-and-y-axes-values?s_tid=srchtitle
    % Loop through each scatter plot
    for i = 1:numel(lineHandles)
        x = lineHandles(i).XData;                         %x data 
        y = lineHandles(i).YData;                         %y data
        xlim(lineHandles(i).Parent, [min(x), max(x)]);    % set x limit to range of x data
        ylim(lineHandles(i).Parent, [min(y), max(y)]);    % set y limit to range of y data
        
        % To convince yourself that the axis scales are still the same within rows/cols,
        % include these two lines of code that will display tick marks.
        %lineHandles(i).Parent.Position(3:4) = lineHandles(i).Parent.Position(3:4) * .8; 
        %set(lineHandles(i).Parent, 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    end
    % now take care of the x axis limits of the histogram plots
    histHandles = h(strcmp(get(h, 'type'), 'histogram'));     %handles to all hist plots
    % loop through hist plots
    for j = 1:numel(histHandles)
        x = histHandles(j).BinEdges;                         %bin edges
        xlim(histHandles(j).Parent, [min(x), max(x)]);       %set x limits
    end
    %%%
    
    title('Correlation of xi'); 
    exportgraphics(gcf, figname); 

end