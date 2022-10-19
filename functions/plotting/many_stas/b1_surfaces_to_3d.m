clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

version_surf = 3; 
lolim = [-94.2132, -62.7473; -86, -72; -88, -76]; 
lalim = [ 51.3480,  18.8802;  30,  44;  36,  32]; 
% i_xsect = 1; 
n_contour = 30; 

depths = [5, 15, 20, 25, 30, 35, 40, ...
        45, 50, 55, 60, 70, 80, 100, 130, 160, 200, 250, 300]; % Try loading these depths. Probably need to type manually for now, but could save as a .mat file in future. 
parms_other = ["zsed", "zmoh"]; 

sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 

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

figure(1); clf; hold on; 
set(gcf, 'pos', [1053 564 767*2 329*ceil(.5*size(lolim,1))], 'color', 'white'); 
tiledlayout(ceil(.5 * size(lolim, 1 )), 2,'TileSpacing', 'Compact')

for i_xsect = 1:size(lolim, 1); 


Q1 = [lalim(i_xsect, 1), lolim(i_xsect, 1)];
Q2 = [lalim(i_xsect, 2), lolim(i_xsect, 2)]; 
[profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));

nxy = 100; 
gcarc = linspace(0, profd, nxy)'; 

[lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, profaz); 


%% Prepare a 2-D section to interpolate into. 


% Simplest section
% lonsect = linspace( min(min(longrid))+1, max(max(longrid))-1, nxy ); 
% latsect = linspace( max(max(latgrid))+1, min(min(latgrid))-1, nxy ); 
zsect = linspace( min(depths), max(depths), nxy-1); 

[lonmesh, zmesh ] = ndgrid(lon_surf_line, zsect); 
[latmesh, zmesh2] = ndgrid(lat_surf_line, zsect); 

mterp = griddata(lon3d, lat3d, z3d, mgrid3d, lonmesh, latmesh, zmesh); 




%% Interpolate sediment, zmoh
zmohsect = griddata(longrid, latgrid, zmoh_surf, lon_surf_line, lat_surf_line); 
zsedsect = griddata(longrid, latgrid, zsed_surf, lon_surf_line, lat_surf_line); 


%%


%% Plot
% figure(1); clf; hold on; 
nexttile(); hold on; 
box on; 
set(gca, 'LineWidth', 1.5, 'YDir', 'reverse'); 
xlabel('Lon'); 
ylabel('Depth (km)'); 
title('Vs cross-section')
colorbar(); 
turbo_map = turbo(); turbo_map = turbo_map(end:-1:1,:); colormap(turbo_map); 

% Velocity
[fk, hand] = contourf(lonmesh, zmesh, mterp, n_contour, 'EdgeAlpha', 0.5); 

% Moho
plot(lon_surf_line, zmohsect, 'LineWidth',4); 
plot(lon_surf_line, zsedsect, 'LineWidth',4); 

end