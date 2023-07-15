clc; clear; restoredefaultpath; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
addpath('/Users/brennanbrunsvik/Documents/repositories/Base_code/colormaps/redblue'); 
tectpath1 = "/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/whitmeyer_karlstrom/layers/"; % Path to what is currently whitmeyer and karlstrom dataset. 
f_distance_pt_to_sta = './_distance_pt_to_sta.mat';
f_xsect_positions = './xsect_positions.mat'; 
pt_dist_nan = 100 /(6371 * 2 * pi / 360); % Don't plot if no station within this many km


%%% Important! Decide which version of figure to plot. Main figure or supplement. 
include_non_velocity = true; % True: show things like radial anisotropy. False: Only show velocity. 

if include_non_velocity; 
    fprintf('\nPlotting main figure paper. \n')
else
    fprintf('\nPlotting supp figure. \n')
end
    

% define colors. 
color_front    = [255, 255, 255]./255; % Grenville Front, Appalachian front.
color_thrust   = [125, 125, 125]./255;
color_mag      = [190, 000, 000]./255;
color_grav     = [000, 150, 000]./255;
color_ha       = [255, 000, 010]./255; % Harrisonburg Anomaly
color_rift     = [255, 111, 000]./255;
color_text     = .95*[1,1,1];
% clim_crust = [3.5, 4.1]; %brb2023/05/10 using different values from xsections to horizontal slices. 
% clim_mantle = [4.25, 4.75]; 

% Define any colors. 
clr_tectfiles = struct("grv_frt", color_front, "MCR", color_rift, ...
    "Reelfoot", color_rift, "something_province", color_thrust); % Field names have to precicely correspond to file names from this tectonic dataset. 

%%% Cross section stuff, copy to cross-section code. 
version_surf = 7; 
ll_min_max_map = [-89  -72   32   46]; % Map view
xsect_positions = load(f_xsect_positions); 
lolim = xsect_positions.lolim; 
lalim = xsect_positions.lalim; 


n_contour = 30; 

if include_non_velocity; 
    depths = [25, 95, 145]; % Try loading these depths. Probably need to type manually for now, but could save as a .mat file in future. 
else
    depths = [10:30:250]; 
end

parms_other = ["zsed", "zmoh"]; 

sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
mdls = load(fresults).mdls; % For sta lon and lat so we know where to and not to plot. 
ta_data = readtable('TA_station_locations.txt'); % Use TA stations, to decide which parts of our results are well constrained. 
lat_sta_ta = ta_data.Latitude ; 
lon_sta_ta = ta_data.Longitude; 
%% Make 3d lat and lon grids. 
nz = length(depths); 
mgrid3d = zeros(size(latgrid,1), size(latgrid,2), nz); 
lat3d   = zeros(size(latgrid,1), size(latgrid,2), nz); 
lon3d   = zeros(size(latgrid,1), size(latgrid,2), nz); 
z3d     = zeros(size(latgrid,1), size(latgrid,2), nz); 
for iz = 1:nz; 
    lat3d(:,:,iz) = latgrid; 
    lon3d(:,:,iz) = longrid; 
    z3d(:,:,iz) = depths(iz); 
end

for idep = 1:length(depths); 
dep = depths(idep); 

this_inversion = sprintf('vs%1.0f',dep); 
sfsmat2= load(sprintf('%s/surface_values_V%1.0f', this_inversion, version_surf)); mgrid_out = sfsmat2.mgrid_out; 

mgrid3d(:,:,idep) = mgrid_out; 

end

% Load crust, sed
sfsmat2 = load(sprintf('zmoh/surface_values_V%1.0f', version_surf)); 
zmoh_surf = sfsmat2.mgrid_out; 
sfsmat2 = load(sprintf('zsed/surface_values_V%1.0f', version_surf)); 
zsed_surf = sfsmat2.mgrid_out; 
sfsmat3 = load(sprintf('xicr/surface_values_V%1.0f', version_surf)); 
xi_surf = sfsmat3.mgrid_out; 

% Load topo
[lon_top, lat_top, z_top...
    ] = get_z_etopo1(min(longrid(:))-1, max(longrid(:))+1, ...
             min(latgrid(:))-1, max(latgrid(:))+1, 'plot', false,...
             'ndmesh', true); 


%% Any borders to plot
%brb2023.02.21 These were copied from the 2021 ENAM paper. They were made only roughly, so I need more accurate shape files. 
app_bord = [-73.9819443, -74.619146 , -75.520073 , -76.1132732, -76.7944389,-77.1789631, -77.6294266, -77.9920049, -78.2117582, -78.4383114,-79.1416076, -79.8888598, -80.2331436, -80.7605221, -81.2438831,-81.9140786, -82.4853569, -83.2434644, -84.3284171, -84.7700317,-85.1106908, -85.5944904, -86.0889287, -87.330563; 40.4970924,  40.7306085,  40.9467137,  41.1455697,  41.2695495,41.2365112,  41.0959121,  40.8719878,  40.7056279,  40.3967643,39.5548831,  38.4965935,  38.1259146,  37.7272803,  37.5445773,37.3439591,  37.1603165,  36.8268747,  36.1733569,  35.880149 ,35.4427709,  34.7777158,  34.1799976,  32.9349287]; 
app_bord2 = [-73.9819  -74.1099  -73.8461  -73.3406  -73.0769  -72.4395; 40.4971   40.9965   42.4397   43.8504   44.8247   45.6140];
app_bord = [flipud(app_bord')', app_bord2]; 
gre_bord = [-82.8790832, -83.5164453, -83.6483134, -84.0878735, -84.3516096, -85.714246 , -87.3406184, -88.1098487; 43.9928145,  41.8040781,  40.8969058,  39.6733704,  38.9764925, 36.6155276,  34.9940038,  34.3978449]; 


%%


ll_min_max_map = [-89  -68   26   46]; % Map view

% Figure out distance from each point to nearest TA station
pt_dist = zeros(size(xgrid)); 
pt_dist_TA = zeros(size(xgrid)); 
for ipt = 1:(size(xgrid, 1) * size(xgrid,2)); 
    pt_dist(ipt) = min(distance(latgrid(ipt), longrid(ipt), ...
        mdls.lat, mdls.lon));
    pt_dist_TA(ipt) = min(distance(latgrid(ipt), longrid(ipt), ...
        lat_sta_ta, lon_sta_ta));
end

[C,~] = contour(longrid, latgrid, pt_dist_TA, [pt_dist_nan, pt_dist_nan]); % Distance where start to nan things before plotting. Distance to TA stations, since those are more what controls our deep resolution 
lon_buff = C(1,2:end)'; % First point is junk for some reason. Remove it. 
lat_buff = C(2,2:end)'; 
rm_buff = (lon_buff > -83) & (lon_buff < -80) & (lat_buff > 44) & (lat_buff < 46); % Manually remove some parts of lon_buff and lat_buff that aren't useful
lon_buff = lon_buff(~rm_buff); 
lat_buff = lat_buff(~rm_buff); 
figure(1010); clf; hold on; 
scatter(lon_buff, lat_buff); 
title('Plot buffer region')

pt_dist_km = pt_dist .* 6371 * 2 * pi / 360; 
pt_dist_TA_km = pt_dist_TA .* 6371 * 2 * pi / 360; 
save(f_distance_pt_to_sta, 'pt_dist_km', 'pt_dist', 'pt_dist_TA', 'lon_buff', 'lat_buff'); % Save for loading in plot_xsects file. 

n_fig_notvel = 3; % How many figures are not velocity but things like moho
n_plots = length(depths) + include_non_velocity * n_fig_notvel; % Add 3 if plotting things that arent just velocity


if include_non_velocity; % Main figure spacing
    figure(17); clf; hold on; set(gcf,'pos', [87 856 692 476]); 
    tiledlayout(2, 3, 'TileSpacing', 'tight'); 
else % Figure spacing if only showing velocity
    figure(17); clf; hold on; set(gcf,'pos', [87 856 692 476*3/2]); 
    tiledlayout(3, 3, 'TileSpacing', 'tight'); 
end

for ifig = 1:n_plots; 

nexttile(); hold on; box on; set(gca, 'LineWidth', 1.5);
m_proj('mercator', 'long',[ll_min_max_map(1), ll_min_max_map(2)],...
                   'lat',[ll_min_max_map(3), ll_min_max_map(4)]);

% Pick the right surface.
cnt_num = 100; 

%% Can plot some things that aren't velocity
if (ifig > n_fig_notvel) || ~ include_non_velocity; 
    ind_depth = ifig; 
    ind_depth = ind_depth - include_non_velocity * n_fig_notvel; % If first few plots aren't velocity, then reduce depth index. 
    grd_cont = mgrid3d(:,:,ind_depth); 

    for ipt = find(~isnan(grd_cont))'; % If any of these within X km of a nan point, kill it. 
        % Compare to TA stations, not just the ones we used. 
        dist_buf = distance(latgrid(ipt), longrid(ipt), lat_buff, lon_buff ); % Distance from the buffered edge of stations to each grid point. Output in degrees, to comapare with pt_dist_nan in degrees
        plt_dist_mod = depths(ind_depth) * 360/(2*pi*6371) ; %pt_dist_nan; % Compare distance, in degrees. 
        if any(dist_buf < plt_dist_mod); 
            grd_cont     (ipt) = nan;
            options.vgrid(ipt) = nan; 
        end
    end

    label = string(sprintf('Z=%1.0f km', depths(ind_depth))) + newline + "Vs (km/s)"; 
    cmap = turbo(12); 
    cmap = flip(cmap); 
    colormap(gca, cmap); 
    if ind_depth == 1; 
        clim([3.6, 3.9]); 
    elseif any(ind_depth == [2,3]); 
        clim([4.35, 4.75]); 
    end
elseif (ifig == 1) && include_non_velocity; % If we want non-velocity plots, then put then in spots 1, 2, 3
    grd_cont = zsed_surf;
    label = "Sediment (km)"; 
    colormap(gca, viridis()); 
    clim([0, 1]);
elseif (ifig == 2) && include_non_velocity; 
    grd_cont = zmoh_surf; 
    label = 'Moho (km)'; 
    colormap(gca, viridis(12)); 
    clim([-25, 50]);    
    clim_min = 5*floor(min(zmoh_surf, [], 'all')/5)  ; % Round in 5s
    clim_max = 5*floor(max(zmoh_surf, [], 'all')/5)  ;
    clim([clim_min, clim_max]); 
elseif (ifig == 3) && include_non_velocity; 
    grd_cont = xi_surf; % TEMPORARY 
    label = '\xi=(Vsh/Vsv)^2'; 
    colormap(gca, redblue(13)); 
    clim(1+[-1, 1].*.1); 
% elseif ifig == 4; 
%     grd_cont = xi_surf; % TEMPORARY 
%     label = 'Vp/Vs'; 
% %     colormap(gca, redblue(13)); 
% %     clim(1+[-1, 1].*.1); 
end

% Make a buffer where we are supposed to nan the model 
% buff_lon, buff_lat

% polyout = polybuffer([buff_lon, buff_lat],'lines', 100 * 360 / (2*pi*6371) ); 
% plot(polyout.Vertices(:,1), polyout.Vertices(:,2))

% Remove stuff from grid outside station area. 
grd_cont(pt_dist > pt_dist_nan) = nan; 
options.vgrid(pt_dist > pt_dist_nan) = nan; 

% for ipt = find(~isnan(grd_cont))'; % If any of these within X km of a nan point, kill it. 
%     % Compare to TA stations, not just the ones we used. 
%     dist_buf = distance(latgrid(ipt), longrid(ipt), lat_buff, lon_buff ); % Distance from the buffered edge of stations to each grid point. Output in degrees, to comapare with pt_dist_nan in degrees
%     plt_dist_mod = 0 ; %pt_dist_nan; 
%     if any(dist_buf < plt_dist_mod); 
%         grd_cont     (ipt) = nan;
%         options.vgrid(ipt) = nan; 
%     end
% end

% Colorbar and label
ax = gca(); 
cbar = colorbar('Location', 'south'); 
cbar.Position(3) = ax.Position(3) * .4; 
cbar.Position(1) = (ax.Position(1)+ax.Position(3)*.5) ; 
if ifig == 1; 
    cbar.Position(1) = cbar.Position(1) - 0.01; % For some reason first axis is off. 
end

% Main contours
num_cnt = 100; 
m_contourf(longrid, latgrid, grd_cont, cnt_num,...
    'LineStyle','none'); 

% State lines
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

% Coast lines
cst = m_coast('patch',[1 1 1], 'FaceAlpha', 0); 

% % % % Cross section position, labels
% % % if ifig == 1; 
% % %     xsect_letters = ["A", "B", "C", "D", "E", "F", "G", "H"]; 
% % %     for ixsect = 1:size(lolim,1); 
% % %         Q1 = [lalim(ixsect, 1), lolim(ixsect, 1)];
% % %         Q2 = [lalim(ixsect, 2), lolim(ixsect, 2)]; 
% % %         [profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));
% % %         gcarc = linspace(0, profd, 100)'; 
% % %         d2km = 2 * pi * 6371 / 360; % Yes I know this is 111, but might as well be somewhat precise :) 
% % %         dist_arc = d2km * gcarc; % km 
% % %         [lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, profaz); 
% % %         m_plot(lon_surf_line, lat_surf_line, 'white', 'linewidth', 1); 
% % %         if ifig == 1; 
% % %             line_shift = 7; 
% % %             color_xsectlab = .8.*[1,1,1]; 
% % %             m_text(lon_surf_line(1+line_shift ), lat_surf_line(1+line_shift  ), xsect_letters(ixsect), ...
% % %                 'color', color_text, 'fontweight', 'normal', 'units', 'data', ...
% % %                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
% % %             m_text(lon_surf_line(end-line_shift), lat_surf_line(end-line_shift)  , xsect_letters(ixsect)+"'", ...
% % %                 'color', color_text, 'fontweight', 'normal', 'units', 'data', ...
% % %                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
% % %         end
% % %     end
% % % end



% % Stations
% if ifig == 2; 
%     m_scatter(mdls.lon, mdls.lat, 5, 'k', 'filled'); 
% end

% Tectonic stuff

% Tectonic fronts. 
m_plot(app_bord(1,:), app_bord(2,:), 'linewidth', 1.5, 'color', color_front); 
% m_plot(gre_bord(1,:), gre_bord(2,:), 'linewidth', 3, 'color', color_front); 

% Load tectonic files to plot. 
plt_ftrs = ["grv_frt"]; % Which things to plot. Not showing thrusts because not pertinent. Those they are mislabeled as mislabeled as  "something_province" % , "MCR", "Reelfoot"
plt_ftrs_path = tectpath1 + plt_ftrs + "*.txt";
for iftr = 1:length(plt_ftrs); 
    fplot = ls(plt_ftrs_path(iftr)); 
    fplot = splitlines(fplot(1:end-1)); % Need to remove last line, to prevent nonsense empty extra cell index. 
    this_color = clr_tectfiles.(plt_ftrs(iftr)); 
    for ifplot = 1:length(fplot); 
        xy = load(fplot{ifplot});  
        x = xy(:,1); 
        y = xy(:,2); 
        m_plot(x, y, 'Color', this_color, 'linewidth', 1.5); % TODO color for each feature
    end
end

% Label features. 
if ifig == 1; 
    m_text(-83, 39, 'GF', 'color', color_text, 'fontsize', 12, ...
        'rotation', 75, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold'); 
    m_text(-83, 36, 'AF', 'color', color_text, 'fontsize', 12, ...
        'rotation', 45, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold'); 
end

% Label the plot
m_text(-74, 30.5, label, 'Units','data', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom'); 

% Handle ticks across different subplots. 
xtck = [-88:6:-68]; 
ytck = [28:4:44]; 
xtckstr = xtck; 
ytckstr = ytck; 

% Ticks for 3x2
if include_non_velocity % Main paper
    xtckcel = {[], [], [], xtckstr, xtckstr, xtckstr}; 
    ytckcel = {ytckstr, [], [], ytckstr, [], []}; 
else % Supplemental figure
    xtckcel = {xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr}; % Manually make 3x3. Sorry, future self. 
    ytckcel = {ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr}; 
end

% Grid, coastlines I think, etc. Background colors. 
m_grid('box','fancy','linestyle','none','gridcolor',.5 .*[1,1,1],...
    'backcolor',[.8, .8, .8], 'xticklabel', xtckcel{ifig}, 'yticklabel', ytckcel{ifig},... % 'backcolor',[.3 .75 1]
    'xtick', xtck, 'ytick', ytck);

% Figure subletter labels - do this late in script so letters don't get covered. 
subplot_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']; 
txt= text(0.04, 0.015, subplot_letters(ifig), 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment','bottom', ...
    'FontSize',13, 'Color',[0,0,0]); 

end

fname_base = 'sage_gage/map_view_V%1.0f%s.%s'; 
if include_non_velocity; 
    sup_txt = ''; 
else
    sup_txt = '_supplement'; 
end

exportgraphics(gcf, sprintf(fname_base, version_surf, sup_txt, 'jpeg'), ...
    'Resolution',500); 
savefig(gcf, sprintf(fname_base, version_surf, sup_txt, 'fig')); 

% exportgraphics(gcf, sprintf('sage_gage/map_view_V%1.0f.jpeg', version_surf), ...
%     'Resolution',500); 
% savefig(gcf, sprintf('sage_gage/map_view_V%1.0f.fig', version_surf)); 

% save()