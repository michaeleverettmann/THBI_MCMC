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

desired_chains = 12; 
desired_iter   = 16000; 

addpath('~/MATLAB/m_map');
% addpath('~/Documents/MATLAB/BayesianJointInv/functions');
addpath('~/Documents/UCSB/ENAM/THBI_ENAM/functions'); 
addpath('~/MATLAB/seizmo/cmap'); warning('adding cmap in seismo. Is this breaking split?'); 
addpath('~/MATLAB/borders'); 
addpath('~/Documents/repositories/general_data'); % For topography loading. get_z_etopo1.m
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 


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



% % % % 
% % % % 
% % % % 
% % % % %% load project & station details
% % % % stations = load([infodir 'stations.mat']); 
% % % % stainfo = stations.stainfo; 
% % % % stainfo.overallQ = ones(size(stainfo.slons)); 
% % % % 
% % % % % % % % retrieve comparison models
% % % % % % % semPath = '~/Documents/repositories/data/models_seismic/SEMum2_avg_VS'
% % % % % % % addpath(semPath); 
% % % % % % % a = SEMum2_avgprofiles(0,[semPath '/']);
% % % % 
% % % % mdls = struct('nwk', {}, 'sta', {}, 'lat', [], 'lon', [], 'model', {},...
% % % %     'misfit', {}, 'dir', {});  % Fill out structure with models we have. 
% % % % 
% % % % 
% % % % %% find stas with good fits overall and good fit to Sp data;
% % % % gdstas = zeros(stainfo.nstas,1);
% % % % for is = 1:stainfo.nstas 
% % % % 
% % % %     clear sdtyp
% % % % 
% % % %     %%%
% % % %     sta = stainfo.stas{is}; 
% % % %     nwk = stainfo.nwk {is}; 
% % % %     
% % % %     resdir = sprintf('%s%s_%s_dat%.0f/%s',proj.STAinversions,sta,nwk,generation,STAMP);
% % % %     fdir = [resdir, '/final_model.mat']; 
% % % %     misfitdir = sprintf('%s/final_misfit.mat',resdir); 
% % % %     gdstas(is) = logical(exist(fdir, 'file'          ) ...
% % % %         * logical(exist(misfitdir, 'file')           ) ...
% % % %         * logical(exist([resdir,'/par.mat'], 'file') )     );  
% % % %     if ~gdstas(is); continue; end; % can't load par if it doesn't exist. Just do continue. 
% % % %     par=load([resdir,'/par.mat']); par=par.par; 
% % % %     gdstas(is) = gdstas(is) * ((par.inv.niter==desired_iter) && (par.inv.nchains==desired_chains)); % Don't plot results if we didn't run the inversion with enough chains or iterations. Might have been test runs..
% % % %     if ~gdstas(is); continue; end; % can't load par if it doesn't exist. Just do continue. 
% % % % 
% % % %     model = load(fdir).final_model; 
% % % %     final_misfit = load(misfitdir).final_misfit; 
% % % % 
% % % %     igdsta = sum(gdstas); 
% % % % 
% % % %     %%% Save useful things to .mat structure. 
% % % %     mdls(1).nwk{igdsta,1} = nwk; 
% % % %     mdls.sta{igdsta,1} = sta; 
% % % %     mdls.lat(igdsta,1) = stainfo.slats(is); 
% % % %     mdls.lon(igdsta,1) = stainfo.slons(is); 
% % % %     mdls.model{igdsta,1} = model; 
% % % %     mdls.misfit{igdsta,1} = final_misfit; 
% % % %     mdls.dir{igdsta,1} = resdir; 
% % % %     %%%
% % % % 
% % % % end
% % % % 
% % % % 
% % % % save(fresults,"mdls"); 