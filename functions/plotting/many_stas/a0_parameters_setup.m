% Compile all the different model files into one. 
% Takes a while to run (a minute or two). 

clear;close all
run("../../../a0_STARTUP_BAYES.m"); 

%% Setup
ifsave = true;

paths = getPaths(); 
paths.STAinversions = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_collate/'; % Place where your results are. 
figPath = '~/Documents/UCSB/ENAM/THBI_ENAM/figures/many_stas/';
infodir = '~/Documents/UCSB/ENAM/THBI_ENAM/ENAM/INFO/'; 
out_dir = '~/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/'; 
feature_dir = '~/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/whitmeyer_karlstrom/layers'; 

desired_chains = 12; 
desired_iter   = 16000; 

addpath('~/MATLAB/m_map');
% addpath('~/Documents/MATLAB/BayesianJointInv/functions');
addpath('~/Documents/UCSB/ENAM/THBI_ENAM/functions'); 
addpath('~/MATLAB/seizmo/cmap'); warning('adding cmap in seismo. Is this breaking split?'); 
addpath('~/MATLAB/borders'); 
addpath('~/Documents/repositories/general_data'); % For topography loading. get_z_etopo1.m
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 
addpath('~/Documents/repositories/Base_code/colormaps/redblue'); 



% specify details of this run
generation = 1; % generation of solution and data processing
STAMP = 'standard';

% Quality thresholds for including stations - important!
overallQ_thresh = 1; % 2 is good, 1 is ok, 0 is bad
Sp_Q_thresh = 1; % Sp data quality (same bounds as above)

if ~ exist(figPath, 'dir'); mkdir(figPath); end
proj = load('~/Documents/UCSB/ENAM/THBI_ENAM/ENAM/INFO/proj.mat').proj; 
proj.STAinversions = paths.STAinversions; 

fresults = sprintf('%s/compiled_results_%s.mat',out_dir,STAMP); 

z_vs = [5:5:300]; warning('Make sure z vs is supposed to be 5:5:300 always in pdfs or something. remove this from parameter setup. ')% Temporary. Different depths looking for station median values. Needs to match whats in the pdfs. 
