function ax = geogax(ax)
% function to make axes with state bounds etc.
    if nargin < 1
        figure(656);clf
        ax = gca;
    end
    axes(ax);
    addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord)
        line(lonbord{iplace}',latbord{iplace}','LineWidth',1.2,'color',0*[1 1 1])
    end
    xlim(ax,[-88 -68]);
    ylim(ax,[30,46]);
    view(ax,0,90)
    set(ax,'fontsize',15,'linewidth',2,'box','on','layer','top','Xgrid','off','Ygrid','off')
end