% brb2023/08/12 Script to gather and save MCMC files for publication in Zenodo
% repository. Most of this was taken from b1_plot_paper_maps_w_litho_mld.m


clc; clear; restoredefaultpath; 
cd .. % brb20240614 Moved to misc folder, so adding cd .. in case the starting folder was important. 
file_upload = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/paper_submission/'; % Name of file we will upload. 
path_many_stas = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/'; 
out_dir = '~/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/'; 

version_surf = 7; 
generation = 1; % generation of solution and data processing
STAMP = 'standard';
z_vs = [5:5:300];
fresults = sprintf('%s/compiled_results_%s.mat',out_dir,STAMP); 


cd(path_many_stas); 

depths = z_vs; 
parms_other = ["zsed", "zmoh"]; 

sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
mdls = load(fresults).mdls; % For sta lon and lat so we know where to and not to plot. 
ta_data = readtable('TA_station_locations.txt'); % Use TA stations, to decide which parts of our results are well constrained. 

%% Get distances from each point to nearest station. 
pt_dist_nan = 100 /(6371 * 2 * pi / 360); % Don't plot if no station within this many km

mdls = load(fresults).mdls;
pt_dist = zeros(size(longrid)); 
for ipt = 1:(size(longrid, 1) * size(longrid,2)); 
    pt_dist(ipt) = min(distance(latgrid(ipt), longrid(ipt), ...
        mdls.lat, mdls.lon));
end

to_nan = pt_dist > pt_dist_nan; 

%% Load velocity and nan outside some distance 
lat_sta_ta = ta_data.Latitude ; 
lon_sta_ta = ta_data.Longitude; 
% Make 3d lat and lon grids. 
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
mgrid_out(to_nan) = nan; 
mgrid3d(:,:,idep) = mgrid_out; 

end

%% Load crust, sed
sfsmat2 = load(sprintf('zmoh/surface_values_V%1.0f', version_surf)); 
zmoh_surf = sfsmat2.mgrid_out; 
sfsmat2 = load(sprintf('zsed/surface_values_V%1.0f', version_surf)); 
zsed_surf = sfsmat2.mgrid_out; 
sfsmat3 = load(sprintf('xicr/surface_values_V%1.0f', version_surf)); 
xi_surf = sfsmat3.mgrid_out; 

zmoh_surf(to_nan) = nan; 
zsed_surf(to_nan) = nan; 
xi_surf  (to_nan) = nan; 


%% Zach's files
T_lab = 1150;

load('zach_lithosphere_results/LAB_MLD/fastlith_map.mat');
load(['zach_lithosphere_results/LAB_MLD/LAB_T',num2str(T_lab),'.mat']);
load('zach_lithosphere_results/LAB_MLD/MLD_vgrads.mat');

%% Plots to make sure things look good.
mcmc_zenodo = struct('lon', lon3d, 'lat', lat3d, 'z', z3d, 'vs', mgrid3d, ...
    'z_moho', zmoh_surf, 'z_sediment', zsed_surf, ...
    sprintf('z_lab_%1.0f', T_lab),z_lab_Tiso_smth); 

save([file_upload 'mcmc_zenodo.mat'], 'mcmc_zenodo'); 

figure(); clf; hold on; 
tiledlayout(2,2, 'TileSpacing','compact'); 
nexttile(); contourf(longrid, latgrid, zmoh_surf); colorbar(); title('zmoh')
nexttile(); contourf(longrid, latgrid, zsed_surf); colorbar(); title('zsed')
nexttile(); contourf(longrid, latgrid, z_lab_Tiso_smth); colorbar(); title('lab')
nexttile(); contourf(longrid, latgrid, mgrid_out); colorbar(); title('Velocity(some depth)')