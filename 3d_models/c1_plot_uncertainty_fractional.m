% Figure S16 in Brunsvik et al 2024. Evaluates the statistical significance
% of anomalies throughout the model. 

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
n_plots = 6; 

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

    %% First, load mean values so we can make a background model. 
    if ifig > fig_offset; 
        dat_type = 'vs'; 
    elseif ifig == 1; 
        dat_type = 'moho'; 
    elseif ifig == 2; 
        dat_type = 'xi'; 
    end

    for i_station = 1:n_stations; 
        if strcmp(dat_type, 'vs'); 
            % Gather data from each station
            param_i = 'VSav'; 
            model = mdls.model{i_station}; 
            parm = model.(param_i); 
            
            zoffset = abs(model.Z - depths(ifig-fig_offset)); 
            [zoffset,idepth] = min(zoffset); 
            zactual = model.Z(idepth); 
    
            parm_all(i_station) = parm(idepth,:); 
        elseif strcmp(dat_type, 'moho'); 
            param_i = 'zmohav'; 
            model =  mdls.model{i_station}; 
            parm_all(i_station) = model.(param_i); 
        elseif strcmp(dat_type, 'xi'); 
            param_i = 'xicrav'; 
            model =  mdls.model{i_station}; 
            parm_all(i_station) = model.(param_i); 
        end; 
    end

    %% Now make a very smooth background model. 
    figure(11); clf; hold on; 
    if strcmp(dat_type, 'xi'); 
        parm_simp = ones(size(parm_all)); 
    else 
        parm_simp = nan(size(parm_all)); 
        for ista = 1:length(parm_simp); 
            dist = distance(lat(ista), lon(ista), lat, lon) * 111; % Distance from this station to each other station
            gauss_dev = 500; 
            weights = exp(-((dist).^2)/(2*gauss_dev^2)); 
            parm_simp(ista) = sum(weights .* parm_all) / sum(weights); 
        end
    end
    scatter(lon, lat, 20, parm_simp, 'filled'); 

    %% Now, figure out what percentage of models (from the ensemble) had a value greater or less than this background. 
    frac_above = nan(size(parm_all)); 
    for ista = 1:length(parm_simp); 
        f_pst = sprintf('%s%s_%s_dat1/standard/posterior.mat',paths.STAinversions, mdls.sta{ista}, mdls.nwk{ista}); 
        pst = load(f_pst).posterior; 

        % Repeat the logic above to extract the posterior ensemble for this parameter.
        if strcmp(dat_type, 'vs'); 
            param_i = 'VSav'; 
            zoffset = abs(pst.zatdep - depths(ifig-fig_offset)); 
            [zoffset,idepth] = min(zoffset); 
            zactual = pst.zatdep(idepth); 
            parm = pst.VSmantle(:,idepth); 
        elseif strcmp(dat_type, 'moho'); 
            parm = pst.zmoh; 
        elseif strcmp(dat_type, 'xi'); 
            parm = pst.xicrust; 
        end; 
        frac_above = sum((parm >= parm_simp(i_station)))/length(parm);        
        parm_all(ista) = frac_above; % For utilizing old script that used parm_all 
    end

    %%
    ax = each_ax(ifig); 
    axes(ax); 
    hold on; box on; set(gca, 'LineWidth', 1.5);
    m_proj('mercator', 'long',[ll_min_max_map(1), ll_min_max_map(2)],...
                       'lat',[ll_min_max_map(3), ll_min_max_map(4)]);

    % For the value to be 1 standard deviation above/below mean, we should have a ratio of 0.841. 

    cutoff = 0.841; 
    toplt = (parm_all > cutoff) | (parm_all < (1-cutoff)); 
    m_scatter(lon( toplt), lat( toplt), 25, parm_all( toplt), 'filled'); 
    m_scatter(lon(~toplt), lat(~toplt), 5 , parm_all(~toplt), 'filled'); 
    clim([0, 1]); 
    colormap(flipud(turbo)); 


    % State lines
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord); 
        m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
    end
    
    % Coast lines
    cst = m_coast('patch',[1 1 1], 'FaceAlpha', 0); 

    % Title 
    m_text(-73.7, 32.8, titles(ifig), 'Units','data', 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom'); 
    
    % Handle ticks across different subplots. 
    xtck = [-88:6:-68]; 
    ytck = [28:4:44]; 
    xtckstr = xtck; 
    ytckstr = ytck; 
    xtckcel = {[], [], [], xtckstr, xtckstr, xtckstr, [],[],[]}; 
    ytckcel = {ytckstr, [], [], ytckstr, [], [], ytckstr, [],[]}; 
    m_grid('box','fancy','linestyle','none','gridcolor',.8 .*[1,1,1],...
        'backcolor',.9.*[1,1,1], 'xticklabel', xtckcel{ifig}, 'yticklabel', ytckcel{ifig},... % 'backcolor',[.3 .75 1]
        'xtick', xtck, 'ytick', ytck);

    % Colorbar and label
    ax = gca(); 
    cbar = colorbar('Location', 'south'); 
    cbar.Position(3) = ax.Position(3) * .4; 
    cbar.Position(1) = (ax.Position(1)+ax.Position(3)*.55) ; 
    cbar.Position(4) = cbar.Position(4) * 2.4; 
    cbar.Position(2) = cbar.Position(2) * 1.01; 

    % Change tick labels to 1, 2, sigma.  
    set(cbar, 'Ticks', [1-0.841, 0.5, 0.841]); 
    set(cbar, 'TickLabels', {'-1\sigma', '0', '1\sigma'}); 

    % Make a histogram of the sigma values. It will overlap the colorbar. 
    ax_temp = axes('Position',cbar.Position); hold on; cla; 
    [pdf_y,pdf_x] = ksdensity(parm_all); 
    plot(pdf_x, pdf_y, 'k', 'LineWidth',3); 
    set(ax_temp, 'Color' , 'None', ...
        'XLim', cbar.Limits, ...
        'XTickLabel', [], 'YTickLabel', [] ); 
    ylabel('P(m>m_0)', 'FontSize',8); 

    disp('One plot done')

end

exportgraphics(gcf, sprintf('sage_gage/map_view_p_exceed_background.jpeg'), ...
    'Resolution',500); 