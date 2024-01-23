clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
addpath('/Users/brennanbrunsvik/Documents/repositories/Base_code/colormaps/redblue'); 
tectpath1 = "/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/whitmeyer_karlstrom/layers/"; % Path to what is currently whitmeyer and karlstrom dataset. 
f_distance_pt_to_sta = './_distance_pt_to_sta.mat';
f_xsect_positions = './xsect_positions.mat'; 

% DEBUG OPTION
downsamp_model = 1; % Downsample the original model by this factor. Only use this for debuging to make griddata faster.  
if downsamp_model ~= 1; 
    warning('DEBUG PLOT. Turn downsamp_model == 1 for your final plot')
end 

% Get lab/mld results. Zach intially sent this in mid/late july. 
lab_mld_fold = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/zach_lithosphere_results/'; % Where Zach sent mld lab results in mid/late july 2023
T_lab = 1150;
lab_file = [lab_mld_fold,'LAB_MLD/LAB_T',num2str(T_lab),'.mat']; 
mld_file = [lab_mld_fold 'LAB_MLD/MLD_vgrads.mat']; 
lab_in = load(lab_file); 
mld_in = load(mld_file); 

pt_dist_nan = 100 /(6371 * 2 * pi / 360); % Don't plot if no station within this many km

% define colors. 
color_front    = [255, 255, 255]./255; % Grenville Front, Appalachian front.
color_thrust   = [125, 125, 125]./255;
color_mag      = [190, 000, 000]./255;
color_grav     = [000, 150, 000]./255;
color_ha       = [255, 000, 010]./255; % Harrisonburg Anomaly
color_rift     = [255, 111, 000]./255;
clim_crust = [3.5, 4.1]; 
clim_mantle = [4.25, 4.75]; 
version_surf = 7; 
max_z = 250; 

% Define any colors. 
clr_tectfiles = struct("grv_frt", color_front, "MCR", color_rift, ...
    "Reelfoot", color_rift, "something_province", color_thrust); % Field names have to precicely correspond to file names from this tectonic dataset. 

%%% Cross section stuff, copy to cross-section code. 
version_surf = 7; 
ll_min_max_map = [-89  -68   32   46]; % Map view
% ll_min_max_map = [-89  -72   32   46]; % Map view
lobase = [-87, -76]; 
labase = [ 43,  35]; 
vdx = 1; 
vdy = 0.8;
dshifts = [-4, -2, 0, 1.5]'; 
lolim = lobase + vdx * dshifts; 
lalim = labase + vdy * dshifts; 

lolim(4,:) = [-83, -70]; 
lalim(4,:) = [42, 43]; 

% Add EW slice
% lolim = [lolim; -88,-74]; 
% lalim = [lalim; 40,40]
% lolim = [lolim; -88,-70]; 
% lalim = [lalim; 40,40]
%%% End cross-section stuff. 

% How far from stations will we plot model? in km. Maybe not in use in this
% script yet? 
min_sta_dist_plt = 50; % Half-distance between TA stations? 

n_contour = 30; 

% depths = [5, 15, 20, 25, 30, 35, 40, ...
%         45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 145, 170, 210, 250, 300]; % Try loading these depths. Probably need to type manually for now, but could save as a .mat file in future. 
depths = [5:5:300]; 

depth_plot = 80; 
idepth_plot = find(depths==depth_plot); 
parms_other = ["zsed", "zmoh"]; 

sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
mdls = load(fresults).mdls; % For sta lon and lat so we know where to and not to plot. 
ta_data = readtable('TA_station_locations.txt'); % Use TA stations, to decide which parts of our results are well constrained. 
lat_sta_ta = ta_data.Latitude ; 
lon_sta_ta = ta_data.Longitude; 

f_3d_model = [pwd() , '/surface_colated_b1_V' , num2str(version_surf) , '.mat']; 

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
sfsmat2= load(sprintf('%s/surface_values_V%1.0f', this_inversion, version_surf)); 
mgrid_out = sfsmat2.mgrid_out; 
% mgrid_out = mgrid_out-mean(mgrid_out(:))+4.3; % SUper rough way to see what would happen if we did dvs
mgrid3d(:,:,idep) = mgrid_out; 

end

% Load crust, sed
sfsmat2 = load(sprintf('zmoh/surface_values_V%1.0f', version_surf)); 
zmoh_surf = sfsmat2.mgrid_out; 
sfsmat2 = load(sprintf('zsed/surface_values_V%1.0f', version_surf)); 
zsed_surf = sfsmat2.mgrid_out; 

% Load topo
[lon_top, lat_top, z_top...
    ] = get_z_etopo1(min(longrid(:))-1, max(longrid(:))+1, ...
             min(latgrid(:))-1, max(latgrid(:))+1, 'plot', false,...
             'ndmesh', true); 


%% Any borders to plot
%brb2023.02.21 These were copied from the 2021 ENAM paper. They were made only roughly, so I need more accurate shape files. 
app_bord = [-73.9819443, -74.619146 , -75.520073 , -76.1132732, -76.7944389,-77.1789631, -77.6294266, -77.9920049, -78.2117582, -78.4383114,-79.1416076, -79.8888598, -80.2331436, -80.7605221, -81.2438831,-81.9140786, -82.4853569, -83.2434644, -84.3284171, -84.7700317,-85.1106908, -85.5944904, -86.0889287, -87.330563; 40.4970924,  40.7306085,  40.9467137,  41.1455697,  41.2695495,41.2365112,  41.0959121,  40.8719878,  40.7056279,  40.3967643,39.5548831,  38.4965935,  38.1259146,  37.7272803,  37.5445773,37.3439591,  37.1603165,  36.8268747,  36.1733569,  35.880149 ,35.4427709,  34.7777158,  34.1799976,  32.9349287]; 
gre_bord = [-82.8790832, -83.5164453, -83.6483134, -84.0878735, -84.3516096, -85.714246 , -87.3406184, -88.1098487; 43.9928145,  41.8040781,  40.8969058,  39.6733704,  38.9764925, 36.6155276,  34.9940038,  34.3978449]; 

%% Save 3-D info, borders, cross-sections, whatever, to plot with Plotly
save(f_3d_model, 'mgrid3d', 'lat3d', 'lon3d', 'z3d', ...
    'latgrid', 'longrid', 'xgrid', 'ygrid', 'depths', ...
    'zmoh_surf', 'zsed_surf'); 

%%
figure(16); clf; hold on; 
% set(gcf, 'pos', [1053 564 767*2 329*ceil(.5*size(lolim,1))], 'color', 'white'); 
% tiledlayout(ceil(.5 * size(lolim, 1 )), 2,'TileSpacing', 'Compact')
% set(gcf, 'pos', [1053 564 767 260*size(lolim,1)])
set(gcf, 'pos', [493 109 674 963]); % Roughly 8.5 by 11
tiledlayout(size(lolim, 1 ), 1,'TileSpacing', 'Compact')

% % % % Load distance from each point to closest station. From b1_plot_paper_maps.m
pt_dist_struct = load(f_distance_pt_to_sta); 
pt_dist = pt_dist_struct.pt_dist_km; 
lat_buff = pt_dist_struct.lat_buff; 
lon_buff = pt_dist_struct.lon_buff; 
% pt_dist_TA = pt_dist_struct.pt_dist_TA; % Degrees
% % % pt_exclude = pt_dist_km > min_sta_dist_plt; 
% % % mgrid3d(pt_exclude) = nan; 
% % % zmoh_surf(pt_exclude) = nan; 
% % % zsed_surf(pt_exclude) = nan; 
% % % % End naning

% Redefine cross-sections so they don't go out of where we have constraints. These variables get erased. 
for ixsect = 1:size(lolim,1); 
    [profd,profaz] = distance(lalim(ixsect,1),lolim(ixsect,1),lalim(ixsect,2),lolim(ixsect,2));
    gcarc = linspace(0, profd, 500)'; 
    [lat_surf_line, lon_surf_line] = reckon(lalim(ixsect, 1), lolim(ixsect,1), gcarc, profaz);
    section_distances    = nan(size(lat_surf_line)); 
    for ipt = 1:length(lat_surf_line); 
        section_distances   (ipt) = min(...
                distance(lat_surf_line(ipt), lon_surf_line(ipt), ...
                mdls.lat, mdls.lon ),...
            [], 'all')*6371*2*pi/360; % Fast enough.
    end
    keep_sect = section_distances < min_sta_dist_plt; 
    keep_sect = find(keep_sect); 
    lolim(ixsect, 1) = lon_surf_line(keep_sect(1)); 
    lolim(ixsect, 2) = lon_surf_line(keep_sect(end)); 
    lalim(ixsect, 1) = lat_surf_line(keep_sect(1)); 
    lalim(ixsect, 2) = lat_surf_line(keep_sect(end)); 
end
save(f_xsect_positions, 'lolim', 'lalim'); 

nxy = 100; 
lat_surf_line_all = zeros(nxy, size(lolim,1) ); 
lon_surf_line_all = lat_surf_line_all; 

% Organize subplots
bord_width = 3; 
buffplt = 0.05; 

n_sections = size(lolim, 1) ; 
ax_x  = [0.07 ] * ones([1, n_sections]); 
% ax_y  = [0.7 ]; 
% ax_y  = [.2:0.2:.8]
ax_y = 0.1 + .190*[0:length(ax_x)-1]; 
ax_dy = [0.15] * ones([1, n_sections]); 
ax_dx = nan(size(ax_dy)); 

% ax_dx = % Not known until width of section is defined. 
% ax_dx = [0.8]; 

% for i_xsect = 1; 
for i_xsect = 1:n_sections;

% if i_xsect ~= 4; 
%     disp('Skipping all but last cross section')
%     continue
% end

%%% Pre-plot prep
% Conversions from xy and latlon
Q1 = [lalim(i_xsect, 1), lolim(i_xsect, 1)];
Q2 = [lalim(i_xsect, 2), lolim(i_xsect, 2)]; 
[profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));
gcarc = linspace(0, profd, nxy)'; 
d2km = 2 * pi * 6371 / 360; 
dist_arc = d2km * gcarc; % km , 1-d
[lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, profaz); 
lat_surf_line_all(:,i_xsect) = lat_surf_line; 
lon_surf_line_all(:,i_xsect) = lon_surf_line; 

% Prepare a 2-D VS section for interpolation
zsect = linspace( min(depths), max(depths), nxy-1); 
[lonmesh, zmesh ] = ndgrid(lon_surf_line, zsect); 
[latmesh, ~] = ndgrid(lat_surf_line, zsect); 
[gcmesh , ~] = ndgrid(gcarc        , zsect); 

% Convert to new evenly spaced basis. downsamp_model option is provided only to make developing code easier, because griddata can take >1min to run. 
mterp = griddata(lon3d   (1:downsamp_model:end,1:downsamp_model:end,1:downsamp_model:end), ...
    lat3d                (1:downsamp_model:end,1:downsamp_model:end,1:downsamp_model:end), ...
    z3d                  (1:downsamp_model:end,1:downsamp_model:end,1:downsamp_model:end), ...
    mgrid3d              (1:downsamp_model:end,1:downsamp_model:end,1:downsamp_model:end), ...
    lonmesh, latmesh, zmesh); % Interpolated section

% % % % Make a cone to cut off stations far from TA. Find distance to TA stations. 
buffer_distance = zeros(size(lat_surf_line)); 
for ipt = 1:length(lat_surf_line); 
    buffer_distance   (ipt) = min(...
            distance(lat_surf_line(ipt), lon_surf_line(ipt), ...
            lat_buff, lon_buff),...
        [], 'all')*6371*2*pi/360; % Fast enough.
end
buffer_keep = zsect < buffer_distance; 
mterp(~buffer_keep) = nan; 

% Station distance and position along line
slon = mdls.lon; 
slat = mdls.lat;
sta_x_ind = zeros(size(slon)); 
sta_dists = zeros(size(slon)); 
sta_cutoff = 70; 
for ista = 1:length(slon); 
    [stadist] = distance(slat(ista), slon(ista), lat_surf_line, lon_surf_line);
    % stadist = stadist 
    close_point = find(min(stadist) == stadist); 
    stadist = stadist(close_point);
    stadist = stadist * 6371 * 2 * pi / 360; 
    sta_dists(ista) = stadist; 
    sta_x_ind(ista) = close_point; 
end
plt_sta_bool = sta_dists <= sta_cutoff; 
sta_dists = sta_dists(plt_sta_bool); 
sta_x_ind = sta_x_ind(plt_sta_bool); 
sta_x_dist = dist_arc(sta_x_ind); 

dist_sect = gcmesh*d2km; % For plotting along section. 2-d

% Interpolate sediment, zmoh, topography
zmohsect = griddata(longrid, latgrid, zmoh_surf, lon_surf_line, lat_surf_line); 
zsedsect = griddata(longrid, latgrid, zsed_surf, lon_surf_line, lat_surf_line); 
ztopsect = interpn (lon_top, lat_top, z_top    , lon_surf_line, lat_surf_line, 'cubic'); % can use interpn here for speed because topo is on a grid. But the surface I inverted isn't on a lon/lat grid (it's a linearly spaced grid in x and y), so we have to use griddata above. 
top_scale = 1/50; 
ztopsect_scaled = -(ztopsect*top_scale)-10; % Scale topography
fprintf('Topography scaling is x%1.0f.\n', top_scale*1000) % e.g.: scale from 1 to 1000 coincident with m to km, then * top_scale = 1/50 to shrink it 50x from there. 

% Set up axes
% ax_box = nexttile(); % Make axis. Change later. 
% ax_box = subplot(3,1,1); hold on; 
axxlims = [min(dist_arc), max(dist_arc)]; 
axylims = [-50, max_z]; 
ax_dx(i_xsect) =  ax_dy(i_xsect) * (diff(axxlims)/diff(axylims)); 
ax_box = axes('Position',[ax_x(i_xsect), ax_y(i_xsect), ...
    ax_dx(i_xsect), ax_dy(i_xsect)]); 

set(ax_box, 'YDir', 'reverse'); 
ax_mantle = copyobj(ax_box, gcf); hold on; % Displayed in middle order.  
ax_crust = copyobj(ax_box, gcf); hold on; 
% linkaxes([ax_box, ax_mantle, ax_crust]); 

%%% Plot velocity
% TODO multiple colors here later
% Color choices. 
% clim_min = 3.9; 
% clim_max = 4.75; 
step_contour = 0.025; 
v_contours = [0:step_contour:10]+0.5*step_contour; 
% clim([clim_min, clim_max]); % TODO temporary 
% n_colors = (clim_max - clim_min) / step_contour; 
n_colors = 25; 

% Prep color map. 
turbo_map = turbo(n_colors); 
turbo_map = turbo_map(end:-1:1,:); 
colormap(turbo_map);
% rbmap = redblue(n_colors); 
% rbmap = rbmap(end:-1:1,:); 
% colormap(rbmap); 


% Contour of velocity
% [fk, hand] = contourf(dist_sect, zmesh, mterp, v_contours, 'EdgeAlpha', 0.1); 

% Mantle contour. Could do a loop. 
ismantle = zmohsect < zmesh; 
iscrust  = zmohsect > zmesh; 
mterpcrust  = mterp; 
mterpmantle = mterp; 
mterpcrust (~iscrust ) = nan; 
mterpmantle(~ismantle) = nan; 
[fk, hand] = contourf(ax_crust , gcmesh*d2km, zmesh, mterpcrust , v_contours, 'EdgeAlpha', 0.1); 
[fk, hand] = contourf(ax_mantle, gcmesh*d2km, zmesh, mterpmantle, v_contours, 'EdgeAlpha', 0.1); 
set(ax_crust , 'ydir', 'reverse'); % For some reason contourf occasionally flips figure to unreverse. 
set(ax_mantle, 'ydir', 'reverse'); % For some reason contourf occasionally flips figure to unreverse. 


%%% Plot other things
moh_color = [250, 2, 192]/250; 
plot(ax_box, gcarc*d2km, zmohsect, 'color', moh_color, 'LineWidth', 4); % Moho
plot(ax_box, gcarc*d2km, ztopsect_scaled, 'k', 'LineWidth', 3); % Topography

%%% Good spot for LAB, MLD
% mld_in, lab_in
% lab_line = interp2(lon_surf_line, )
labsect = griddata(lab_in.longrid, lab_in.latgrid, lab_in.z_lab_Tiso_smth, ...
    lon_surf_line, lat_surf_line); 
plot(ax_box, gcarc*d2km, labsect, 'color', moh_color, 'LineWidth', 4,...
    'LineStyle','--'); % Moho


mld_info = mld_in.MLD_info; 
zmld_pref = mld_info.zmld_pref; 
mldsect = griddata(mld_info.longrid, mld_info.latgrid, zmld_pref, ...
    lon_surf_line, lat_surf_line, 'linear'); 
% plot(ax_box, gcarc*d2km, mldsect, 'color', moh_color, 'LineWidth', 4,...
%     'LineStyle',':'); % Moho
plot(ax_box, gcarc*d2km, mldsect, 'color', 1.*[1,1,1], 'LineWidth', 4,...
    'LineStyle','-'); % Moho
% plot(ax_box, gcarc*d2km, mldsect, 'color', [255,255,255]./255, 'LineWidth', 4,...
%     'LineStyle','-'); % Moho
% scatter(ax_box, gcarc*d2km, mldsect, 25, '+', 'color', 'k', 'linewidth', 0.75); % Moho
%%%

% Plot stations
sta_y_plt = interp1(dist_arc, ztopsect_scaled, sta_x_dist); 
fck = scatter(ax_box, sta_x_dist, sta_y_plt-10, ...
    'marker', '^', 'CData', [0,0,0], 'SizeData', 20, 'LineWidth', 0.25); 
% 'MarkerFaceColor','k', 

% Tectonic borders
for ibord = 1:2
    if ibord == 1; 
        bord = app_bord; 
        textplt = 'AF'; 
    elseif ibord == 2; 
        bord = gre_bord; 
        textplt = 'GF'; 
    end
    [profd_bord,profaz_bord] = distance(Q1(1),Q1(2),bord(2,:), bord(1,:)); % Distance from origin, and azimuth. 
    bord_dist = d2km * interp1(profaz_bord, profd_bord, profaz, 'spline');  % Interpolate for the borders distance where the azimuth matches this sections azimuth. 
    bord_height = interp1(dist_arc, ztopsect_scaled, bord_dist); % Plot at top of topography
    scatter(ax_box, bord_dist, bord_height, ...
        450, 'k', '|', 'linewidth', bord_width); 
%     text(ax_box, bord_dist+13, bord_height+7, textplt, ...
%         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
%         'FontSize',8); 
    text(ax_box, bord_dist+19, -35, textplt, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'FontSize',15, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3); 

end

% % Continent ocean transition, based on topography
plot_COT = false; 
if plot_COT; 
    is_land = find(ztopsect > 0); 
    last_land = is_land(end);
    fprintf('COT determination - assuming ocean is on the high dist_arc side of section.')
    if last_land < length(lon_surf_line); % Ignore if it's the last index
        ocean_x = dist_arc(last_land); 
        ocean_y = ztopsect_scaled(last_land); 
        scatter(ax_box, ocean_x, ocean_y, ...
            450, 'k', '|', 'linewidth', bord_width); 
        text(ax_box, ocean_x - 20, ocean_y-20, 'COT', 'HorizontalAlignment', 'Right', 'FontSize', 11); 
    end
end

% Section letters. 
section_letter = char(64+n_sections+1-i_xsect); % Text for cross-section name. ith letter of alphabet
t1=text(ax_box, 0.01, 1.14, section_letter    , 'fontsize', 18, 'color', 'k', 'units', 'normalized', 'VerticalAlignment','top'); 
t2=text(ax_box, 0.99, 1.14, section_letter+"'", 'fontsize', 18, 'color', 'k', 'units', 'normalized', 'VerticalAlignment','top', 'HorizontalAlignment','right'); 

% Feature labels, very manual 
if i_xsect == 1; 
    text(ax_box, 270, 115, 'SAA', 'fontsize', 16);
elseif i_xsect == 3;
    text(ax_box, 660, 100, 'CAA', 'fontsize', 16);
elseif i_xsect == 4;
    text(ax_box, 820, 100, 'NAA', 'fontsize', 16);
end

% Exageration
dexag = 25; 
plot(ax_box, [min(dist_arc)+dexag, min(dist_arc)+dexag] , [max_z, max_z-dexag], 'k', 'linewidth', 1.5); 
plot(ax_box, [min(dist_arc), min(dist_arc)+dexag] , [max_z-dexag, max_z-dexag], 'k', 'linewidth', 1.5); 


% % % % Manual things
% % % ylim([-70, 300]); 
% % % text(400, -40, 'Grenville front', 'HorizontalAlignment', 'left', 'FontSize', 11); 
% % % text(830, -40, 'Appalachian front', 'HorizontalAlignment', 'left', 'FontSize', 11); 
% % % text(1270, 40, 'Moho', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 13, 'color', moh_color); 


linkaxes([ax_box, ax_mantle, ax_crust]); % done already, but I guess do it again. 
axes(ax_crust); axes(ax_mantle); axes(ax_box); 

% Some axis set up that only works after plots are made for some reason
ax_box.Visible = 'on' ; ax_mantle.Visible = 'off'; ax_crust.Visible = 'off'; 
ax_box.Color = 'none'; ax_mantle.Color = 'none'; ax_crust.Color = 'none'; 
box(ax_box, 'on'); 
set(ax_box, 'LineWidth', 1.5); 
ylim(ax_box, axylims); 
xlim(ax_box, axxlims);
clim(ax_crust , clim_crust ); 
clim(ax_mantle, clim_mantle); 

% Change topo labels
% ax_box.YTickLabel{1} = []; 
ax_box.YTick = [0:50:10000]; % Steps of 50, starting at 0. Nothing for topography. 



ylabel('Depth (km)'); 
if i_xsect == 1; 
    xlabel('Distance (km)'); 
end
set(gca, 'FontSize', 12);

end


% % % % Colorbars if horizontal
% % % cbar_crust  = colorbar(ax_crust ,'Location', 'south'); 
% % % cbar_mantle = colorbar(ax_mantle,'Location', 'south'); 
% % % cbar_crust.Position      = [buffplt     , .035, .5-2*buffplt, .015]; 
% % % cbar_mantle.Position     = [.5 + buffplt, .035, .5-2*buffplt, .015]; 
% % % cbar_crust.Label.String  = 'Crust Vs' ; 
% % % cbar_mantle.Label.String = 'Mantle Vs'; 
% Colorbars if vertical
cbar_crust  = colorbar(ax_crust ,'Location', 'east'); 
cbar_mantle = colorbar(ax_mantle,'Location', 'east'); 
set(cbar_crust , 'FontSize', 12);
set(cbar_mantle, 'FontSize', 12);

% cbar_crust.Position      = [ax_dx(1) + 0.1  , ax_y(1), .02, ax_dy(1)]; 
% cbar_mantle.Position     = [ax_dx(1) + 0.2, ax_y(1), .02, ax_dy(1)]; 
cbar_crust.Position      = [ax_dx(end) + 0.08  , ax_y(end), .0175, ax_dy(1)]; 
cbar_mantle.Position     = [ax_dx(end) + 0.16 , ax_y(end), .0175, ax_dy(1)]; 
cbar_crust.Label.String  = 'Crust Vs' ; 
cbar_mantle.Label.String = 'Mantle Vs'; 

%% Map view of cross-section locations. 
% axmap = axes('Position', [ax_dx(end)+0.075, (ax_y(end)+ax_dy(end)/6), ax_dy(end)*.8, ax_dy(end)*.8]);
axmap = axes('Position', [ax_dx(1)+0.1, (ax_y(1)-0.025), ax_dy(1)*1.25, ax_dy(1)*1.25]);
cla(axmap); hold on; 
% set(gca, 'FontSize', 12);
m_proj('mercator', 'long',[ll_min_max_map(1)-1, ll_min_max_map(2)+1],...
                   'lat',[ll_min_max_map(3)-2, ll_min_max_map(4)]); 
% box on; 
% set(gca,'Visible', 'off'); 
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',.5*[1 1 1])
end
cst = m_coast('patch',.5*[1 1 1], 'FaceAlpha', 0); 


% Appalachian and Grenville fronts. 
app_bord = [-73.9819443, -74.619146 , -75.520073 , -76.1132732, -76.7944389,-77.1789631, -77.6294266, -77.9920049, -78.2117582, -78.4383114,-79.1416076, -79.8888598, -80.2331436, -80.7605221, -81.2438831,-81.9140786, -82.4853569, -83.2434644, -84.3284171, -84.7700317,-85.1106908, -85.5944904, -86.0889287, -87.330563; 40.4970924,  40.7306085,  40.9467137,  41.1455697,  41.2695495,41.2365112,  41.0959121,  40.8719878,  40.7056279,  40.3967643,39.5548831,  38.4965935,  38.1259146,  37.7272803,  37.5445773,37.3439591,  37.1603165,  36.8268747,  36.1733569,  35.880149 ,35.4427709,  34.7777158,  34.1799976,  32.9349287]; 
app_bord2 = [-73.9819  -74.1099  -73.8461  -73.3406  -73.0769  -72.4395; 40.4971   40.9965   42.4397   43.8504   44.8247   45.6140];
app_bord = [flipud(app_bord')', app_bord2]; 
m_plot(app_bord(1,:), app_bord(2,:), 'Color', 'k', 'linewidth', 2, 'linestyle', ':'); % , 'LineStyle', '-.');

% Load Grenville file and plot 
plt_ftrs_path = tectpath1 + "grv_frt*.txt";
fplot = ls(plt_ftrs_path); 
fplot=fplot(1:end-1)
xy = load(fplot);  
m_plot(xy(:,1), xy(:,2), 'Color', 'k', 'linewidth', 2, 'linestyle', ':'); 

% cross section locations
xsect_letters = ["A", "B", "C", "D"]; 
xsect_letters = flip(xsect_letters); 
% xsect_letters = xsect_letters(size(lolim):-1:1); % Flip letters so we label going from top to bottom
for ixsect = size(lolim,1):-1:1; % 1:size(lolim,1); 
    Q1 = [lalim(ixsect, 1), lolim(ixsect, 1)];
    Q2 = [lalim(ixsect, 2), lolim(ixsect, 2)]; 
    [profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));
    gcarc = linspace(0, profd, 100)'; 
    d2km = 2 * pi * 6371 / 360; % Yes I know this is 111, but might as well be somewhat precise :) 
    dist_arc = d2km * gcarc; % km 
    [lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, profaz); 
    m_plot(lon_surf_line, lat_surf_line, 'k', 'linewidth', 2); 
    m_text(lon_surf_line(1 )-.75, ...
        lat_surf_line(1  )-0, ...
        xsect_letters(ixsect), ...
        'color', 'k', 'fontweight', 'normal', 'units', 'data', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'fontsize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 0.5); 
    m_text(lon_surf_line(end) +1.75, ...
        lat_surf_line(end),...
        xsect_letters(ixsect)+"'", ...
        'color', 'k', 'fontweight', 'normal', 'units', 'data', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'fontsize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 0.5);   
end

% m_grid('box','fancy','linestyle','none','gridcolor',.5 .*[1,1,1],...
%     'backcolor','none', 'xtick', [-85:10:75], 'ytick', [35:10:45]); 
m_grid('box','fancy','linestyle','none','gridcolor',.5 .*[1,1,1],...
    'backcolor','none', 'xtick', [-85:10:75], 'ytick', [35:10:45], ...
    'fontsize', 12, 'ylabeldir', 'middle', 'YaxisLocation', 'right'); 

%%

mkdir('xsections'); 
set(gcf, 'Renderer', 'painters'); 
exportgraphics(gcf, sprintf('sage_gage/xsections_V%1.0f.jpeg', version_surf), ...
    'Resolution', 300); 
% savefig(gcf, sprintf('sage_gage/xsections_V%1.0f.fig', version_surf)); 




% %% Plot of cross-section positions
% % ll_min_max_map
% a3_2_plot_surface_simple(llminmax, 'stalon', mdls.lon, 'stalat', mdls.lat, ...
%     'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', mgrid3d(:, :, idepth_plot),...
%     'sectlon', lon_surf_line_all, 'sectlat', lat_surf_line_all); 
% text(0.75, 0.1, sprintf('Z=%1.0f km',depth_plot), 'Units', 'normalized'); 
% set(gcf, 'pos',[295 476 426 321]); 
% exportgraphics(gcf, sprintf('sage_gage/xsections_map_V%1.0f.pdf', version_surf)); 
% fprintf('\n'); 