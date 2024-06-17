run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

% Parameters to set. Taken initially from b1_plot_paper_maps_nat_geo.m
ll_min_max_map = [-89  -72   32   46]; % Map view
ll_min_max_map = [-89  -68   26   46]; % Map view

depths = [25, 60, 95, 145]; 
parms_other = [];%"zsed", "zmoh"]; 
titles = ["Moho", "Xi", "Z="+string(depths)]; 


mdls = load(fresults).mdls; % Get station positions and uncertainties. 

lat = mdls.lat; 
lon = mdls.lon; 

param_plot = [string(depths), parms_other]; 
n_plots = 6; % length(param_plot); 
% i_param = 1; 

%%
figure(1); clf; hold on; set(gcf,'pos', [87 856 692 476*3/2]); 
t = tiledlayout(3, 3, 'TileSpacing', 'tight');

% Make an array of the axes in tiled layout. Then we can call axes(each_ax(ifig)). Calling nexttile in the loop is not compatible with making the sub axis that overlaps the colorbar. 
each_ax = []; 
for ifig = 1:n_plots; 
    each_ax = [each_ax, nexttile(ifig)]; 
end

for ifig = 1:n_plots; 

    n_stations = size(mdls.lat,1); 
    parm_all = nan(n_stations, 1); 

    fig_offset = 2; 
    for i_station = 1:n_stations; 
        if ifig > fig_offset; 
            % Gather data from each station
            param_i = 'VSsig2'; 
            model = mdls.model{i_station}; 
            parm = model.(param_i); 
            
            zoffset = abs(model.Z - depths(ifig-fig_offset)); 
            [zoffset,idepth] = min(zoffset); 
            zactual = model.Z(idepth); 
    
            parm = parm(idepth,:); 
            parm = (parm(2) - parm(1))/2; % NOTE dividing by 2, because I think that we have: (mu+sigma1)-(mu-sigma1)=2*sigma1, 
            parm_all(i_station) = parm; 

        elseif ifig == 1; 
            param_i = 'zmohsig'; 
            model =  mdls.model{i_station}; 
            parm = model.(param_i); 
            parm = parm*2; % Multiply by 2, so it is like sigma2 instead of sigma 1. We do not have sigma2 value saved. 
            parm_all(i_station) = parm; 
        elseif ifig == 2; 
            param_i = 'xicrsig2'; 
            model =  mdls.model{i_station}; 
            parm = model.(param_i); 
            parm = (parm(2) - parm(1))/2; % NOTE dividing by 2, because I think that we have: (mu+sigma1)-(mu-sigma1)=2*sigma1,   
            parm_all(i_station) = parm;
        end; 
    end

    ax = each_ax(ifig); 
    axes(ax); 
    hold on; box on; set(gca, 'LineWidth', 1.5);
    m_proj('mercator', 'long',[ll_min_max_map(1), ll_min_max_map(2)],...
                       'lat',[ll_min_max_map(3), ll_min_max_map(4)]);
    
    m_scatter(lon, lat, 20, parm_all, 'filled'); 

    % State lines
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord); 
        m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
    end
    
    % Coast lines
    cst = m_coast('patch',[1 1 1], 'FaceAlpha', 0); 

    % Title 
    m_text(-74, 32.5, titles(ifig), 'Units','data', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom'); 
    
    % Handle ticks across different subplots. 
    xtck = [-88:6:-68]; 
    ytck = [28:4:44]; 
    xtckstr = xtck; 
    ytckstr = ytck; 
    xtckcel = {[], [], [], xtckstr, xtckstr, xtckstr, [],[],[]}; 
    ytckcel = {ytckstr, [], [], ytckstr, [], [], ytckstr, [],[]}; 
    % Grid, coastlines I think, etc. Background colors. 
    m_grid('box','fancy','linestyle','none','gridcolor',.8 .*[1,1,1],...
        'backcolor',.9.*[1,1,1], 'xticklabel', xtckcel{ifig}, 'yticklabel', ytckcel{ifig},... % 'backcolor',[.3 .75 1]
        'xtick', xtck, 'ytick', ytck);

    % Colorbar and label
    ax = gca(); 
    cbar = colorbar('Location', 'south'); 
%     cbar.Position(3) = ax.Position(3) * .47; 
%     cbar.Position(1) = (ax.Position(1)+ax.Position(3)*.47) ; 
%     cbar.Position(4) = cbar.Position(4) * 2.4; 
    cbar.Position(3) = ax.Position(3) * .4; 
    cbar.Position(1) = (ax.Position(1)+ax.Position(3)*.55) ; 
    cbar.Position(4) = cbar.Position(4) * 2.4; 
    colormap('turbo'); 

    % Make a histogram of the sigma values. It will overlap the colorbar. 
    ax_temp = axes('Position',cbar.Position); hold on; cla; 
    [pdf_y,pdf_x] = ksdensity(parm_all); 
    plot(pdf_x, pdf_y, 'k', 'LineWidth',3); 
    set(ax_temp, 'Color' , 'None', ...
        'XLim', cbar.Limits, ...
        'XTickLabel', [], 'YTickLabel', [] ); 
    ylabel('P(m)'); 
    disp('One plot done')

end


exportgraphics(gcf, sprintf('sage_gage/map_view_sigma.jpeg'), ...
    'Resolution',500); 
% savefig(gcf, sprintf(fname_base, version_surf, sup_txt, 'fig')); 