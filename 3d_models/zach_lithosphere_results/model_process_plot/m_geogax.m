function ax = m_geogax(ax)
% function to make axes with state bounds etc., for m_map plot
    if nargin < 1
        figure(656);clf
        ax = gca;
    end
    axes(ax);
    addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
    addpath('/Users/zeilon/Dropbox/MATLAB/lib/m_map/')
    ll_min_max_map = [ -89   -67    30    46 ];
    m_proj('mercator', 'long',[ll_min_max_map(1), ll_min_max_map(2)],...
                   'lat',[ll_min_max_map(3), ll_min_max_map(4)]); 
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord)
        m_line(lonbord{iplace}',latbord{iplace}','LineWidth',1.2,'color',0*[1 1 1])
    end
%     xlim(ax,[-89 -69]);
%     ylim(ax,[26,46]);
%     view(ax,0,90)
    set(ax,'fontsize',15,'linewidth',2,'box','on','layer','top','Xgrid','off','Ygrid','off')
    m_grid
end