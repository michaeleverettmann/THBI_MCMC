clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

version_surf = 7; 
% lolim = [-87, -76; -86, -68; -88, -78; -87, -80.5]; 
% lalim = [ 43,  35;  30,  47;  36,  33;  38,  25  ]; 
lobase = [-87, -76]; 
labase = [ 43,  35]; 
vdx = 1; 
vdy = 0.8;
% vshift = [vdx, vdy]; 
dshifts = flip([-7, -5, -3, -1, 0, 1, 5, 7]'); 
lolim = lobase + vdx * dshifts; 
lalim = labase + vdy * dshifts; 


% i_xsect = 1; 
n_contour = 30; 

depths = [5, 15, 20, 25, 30, 35, 40, ...
        45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 145, 170, 210, 250, 300]; % Try loading these depths. Probably need to type manually for now, but could save as a .mat file in future. 
parms_other = ["zsed", "zmoh"]; 

sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
mdls = load(fresults).mdls; % For sta lon and lat so we know where to and not to plot. 


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

%%

figure(16); clf; hold on; 
% set(gcf, 'pos', [1053 564 767*2 329*ceil(.5*size(lolim,1))], 'color', 'white'); 
% tiledlayout(ceil(.5 * size(lolim, 1 )), 2,'TileSpacing', 'Compact')
set(gcf, 'pos', [1053 564 767 329*size(lolim,1)])
tiledlayout(size(lolim, 1 ), 1,'TileSpacing', 'Compact')

nxy = 100; 
lat_surf_line_all = zeros(nxy, size(lolim,1) ); 
lon_surf_line_all = lat_surf_line_all; 

% for i_xsect = 1; 
for i_xsect = 1:size(lolim, 1); 



Q1 = [lalim(i_xsect, 1), lolim(i_xsect, 1)];
Q2 = [lalim(i_xsect, 2), lolim(i_xsect, 2)]; 
[profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));

gcarc = linspace(0, profd, nxy)'; 

[lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, profaz); 

lat_surf_line_all(:,i_xsect) = lat_surf_line; 
lon_surf_line_all(:,i_xsect) = lon_surf_line; 


%% Prepare a 2-D section to interpolate into. 


% Simplest section
% lonsect = linspace( min(min(longrid))+1, max(max(longrid))-1, nxy ); 
% latsect = linspace( max(max(latgrid))+1, min(min(latgrid))-1, nxy ); 
zsect = linspace( min(depths), max(depths), nxy-1); 

[lonmesh, zmesh ] = ndgrid(lon_surf_line, zsect); 
[latmesh, zmesh2] = ndgrid(lat_surf_line, zsect); 
[gcmesh , zmesh3] = ndgrid(gcarc        , zsect); 

mterp = griddata(lon3d, lat3d, z3d, mgrid3d, lonmesh, latmesh, zmesh); 




%% Interpolate sediment, zmoh
zmohsect = griddata(longrid, latgrid, zmoh_surf, lon_surf_line, lat_surf_line); 
zsedsect = griddata(longrid, latgrid, zsed_surf, lon_surf_line, lat_surf_line); 


%%


%% Plot
% figure(1); clf; hold on; 
d2km = 2 * pi * 6371 / 360; % Yes I know this is 111, but might as well be somewhat precise :) 
nexttile(); hold on; 
box on; 
set(gca, 'LineWidth', 1.5, 'YDir', 'reverse'); 
xlabel('Distance (km)'); 
ylabel('Depth (km)'); 
title('Vs cross-section')
colorbar(); 
turbo_map = turbo(); turbo_map = turbo_map(end:-1:1,:); colormap(turbo_map); 
clim([3.5, 4.8]); % TODO temporary 

% Velocity
[fk, hand] = contourf(gcmesh*d2km, zmesh, mterp, n_contour, 'EdgeAlpha', 0.5); 

% Moho
plot(gcarc*d2km, zmohsect, 'k', 'LineWidth', 5); 
plot(gcarc*d2km, zsedsect, 'k', 'LineWidth', 5); 

section_letter = char(64+i_xsect); % Text for cross-section name. ith letter of alphabet
t1=text(0.01, 1.12, section_letter    , 'fontsize', 20, 'color', 'r', 'units', 'normalized', 'VerticalAlignment','top'); 
t2=text(0.99, 1.12, section_letter+"'", 'fontsize', 20, 'color', 'r', 'units', 'normalized', 'VerticalAlignment','top', 'HorizontalAlignment','right'); 



end

linkaxes(gcf().Children.Children); 


mkdir('xsections'); 
exportgraphics(gcf, sprintf('xsections/xsections_V%1.0f.pdf', version_surf), 'Resolution', 300); 


%% Plot of cross-section positions
a3_2_plot_surface_simple(llminmax, 'stalon', mdls.lon, 'stalat', mdls.lat, ...
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', mgrid3d(:, :, end-7),...
    'sectlon', lon_surf_line_all, 'sectlat', lat_surf_line_all); 
exportgraphics(gcf, sprintf('xsections/xsections_map_V%1.0f.pdf', version_surf)); 