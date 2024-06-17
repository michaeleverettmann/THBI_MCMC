%% Pseudo cross section for all sites with inversions along section
clear;close all

run("../../a0_STARTUP_BAYES.m"); 

%% Setup
paths = getPaths(); 
proj = load('~/Documents/UCSB/ENAM/THBI_ENAM/ENAM/INFO/proj.mat'); 
proj = proj.proj; 
paths.STAinversions = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_collate/'; % Place where your results are. 
proj.STAinversions = paths.STAinversions; 
sta_list_path = '~/Documents/UCSB/ENAM/THBI_ENAM/ENAM/batch/staList_all.txt'; 

figPath = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/figures/xsect/'; warning('Export path not relative path')
    
addpath('~/MATLAB/m_map');
addpath('~/Documents/UCSB/ENAM/THBI_ENAM/functions'); 
addpath('~/MATLAB/seizmo/cmap'); warning('adding cmap in seismo. Is this breaking split?'); 
addpath('~/MATLAB/borders'); 
addpath('/Users/brennanbrunsvik/Documents/repositories/general_data'); % For topography loading. get_z_etopo1.m

% specify details of this run
generation = 1; % generation of solution and data processing
STAMP = 'standard';

% Quality thresholds for including stations - important!
overallQ_thresh = 1; % 2 is good, 1 is ok, 0 is bad
Sp_Q_thresh = 1; % Sp data quality (same bounds as above)

ifsave = true;


%% lon/lat limits on stations to include:
lolim = [-90, -68]; 
lalim = [26, 45];

ofile1 = [figPath 'srough_sta_map_',STAMP];

infodir = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/INFO/'; 
stations = load([infodir 'stations.mat']); 
stainfo = stations.stainfo; 
stainfo.overallQ = ones(size(stainfo.slons)); 

%% Only look at stas we want to run inversion on
sta_list_all = string(table2cell(readtable(...
        sta_list_path, 'ReadVariableNames', 0))); 
sta_list_all = sta_list_all(:,1) + ' ' + sta_list_all(:,2); % Combine net and sta
sta_list_all = sta_list_all(1:end-1,:); % Remove non station final line

stainfo.netsta = string(stainfo.nwk) + ' ' + string(stainfo.stas); 

gdstas = zeros(length(stainfo.slons),1); 
for igdstas = 1:length(stainfo.slons); 
    if any(stainfo.netsta(igdstas) == sta_list_all); 
        gdstas(igdstas) = true; 
    else; 
        gdstas(igdstas) = false; 
    end
end

gdstas = logical(gdstas); 

%% Making my own map brb2022.03.08
mapFigNum = 1001; 
figure(mapFigNum); clf; hold on; set(gcf, 'color', 'white', 'pos', [-1152 439 378 369]); 
m_proj('lambert','long',lolim + [-2 2],'lat',lalim + [-2 2]);
m_coast('patch',[1 .85 .7]);

[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);

staScat = m_plot(stainfo.slons(gdstas),stainfo.slats(gdstas),...
    '^','linewidth',0.01, 'markerSize', 14, ...
    'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0]); 