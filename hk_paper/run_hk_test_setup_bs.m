% Obnoxious code I don't want to look at. Gets mcmc folder stuff ready, to
% do forward modelling of receiver functions. Also some folder management. 

% close all 
% clear all % Cannot clear all, because then we allways reset network_manual and station_manual below. 
%% Setup
run('../a0_STARTUP_BAYES.m')

addpath('/Users/brennanbrunsvik/Documents/repositories/hk_anis'); %brb TODO need to move this. 

fig_path = sprintf('%s/../figures/hk_paper', pwd()) ; 
fhand_figname = @(zmoh, k, thisfig, frmt)sprintf(...
    '%s/z%3.0f_k%3.0f_%s.%s',fig_path, zmoh*10, k*100, thisfig, frmt); % Convenient function to make figure names. Get rid of decimals. 

proj = struct('name', 'SYNTHETICS'); % bb2021.08.04 changed from EXAMPLE because I don't have the example data files. %,'EXAMPLE');
paths = getPaths(); 
proj.STAinversions = paths.STAinversions; 
proj.dir = [fileparts(mfilename('fullpath')),'/'];
proj.STAinversions = paths.STAinversions; ; % [proj.dir,'inversion_results/'];
save([proj.dir,'project_details.mat'],'proj')

wd = pwd; addpath(wd);
cd(proj.dir);

%% specify details of this run
generation = 0; % generation of solution and data processing
gc = '';
BWclust = '';
STAMP = 'hk_paper_1';
onesta = '';

%% put parameters in place 
global run_params

% you can obviously adapt these (likely in some loop over stations etc.) to be appropriate for your dataset
run_params.projname = proj.name; % from above
run_params.gc = gc; % great circle distance of body wave data, if relevant
run_params.BWclust = BWclust; % cluster of BW data, if relevant
run_params.datN = generation; % processing iteration, if relevant
run_params.STAMP = STAMP; % NEED - some identifier for this inversion run
run_params.overwrite = 1; % do you want to overwrite previous results?
% % % if ~ (exist('network_manual', 'var') && exist('station_manual', 'var')) ; 
network_manual = 'testnwk'; 
station_manual = 'simple_layers_1'; %
fprintf('\nReseting to %s.%s\n',network_manual,station_manual)
% % % end
run_params.sta = station_manual; % name of station
run_params.nwk = network_manual; % name of network


global run_params
paths = getPaths(); 

projname = run_params.projname;
sta = run_params.sta;
nwk = run_params.nwk;
gc = run_params.gc;
BWclust = run_params.BWclust;
datN = run_params.datN;
STAMP = run_params.STAMP;
overwrite = run_params.overwrite;
global projdir TRUEmodel
projdir = [paths.THBIpath,'/',projname,'/'];
cd(projdir);
run([paths.THBIpath,'/a0_STARTUP_BAYES']);
load('project_details'); %TODO_STATION_NETWORK bb2021.11.12
addpath([proj.dir,'matguts/']);

%% PARMS
run parms/bayes_inv_parms
[par, inv] = update_bayes_inv_parms(par, STAMP); % Modify this function to make different tests. 


if strcmp(projname,'SYNTHETICS')
    par.stadeets = struct('sta',sta','nwk',nwk'); 
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise
end

if strcmp(projname,'SYNTHETICS') || strcmp(projname,'LAB_tests')
    par.synth.noise_sta_deets = struct('datadir',['/Volumes/data/THBI/US/STAsinv/',noisesta,'_dat20/'],...
                         'sta',noisesta,'nwk',noisenwk,'gc',noisegcarcs,'noiseup',noiseup,'noiseshape',noiseshape);
end

par.inv.BWclust = BWclust;
ifsavedat = false;

%% get saving things ready
par.proj = proj;
avardir = sprintf('%s%s_%s_dat%.0f/',par.proj.STAinversions,sta,nwk,datN);
resdir = [avardir,STAMP];
if ~exist(resdir,'dir'), try mkdir(resdir); catch, error('Looks like no path to output directory - is it mounted?'); end, end

par.data = struct('stadeets',struct('sta',sta,'nwk',nwk,'Latitude',[],'Longitude',[]),...
                  'gc',gc,'datN',datN,'avardir',avardir);

par.res.STAMP = STAMP;
par.res.resdir= resdir;
par.res = orderfields(par.res,{'STAMP','resdir','zatdep'});

%% Get some directories ready. 
% Switch to execution folder, to make synthetic data. 
prev_dir = pwd(); 
cd(paths.ramDrive); % Execute everything from a folder in ram for major speedup. 
mkdir([nwk '_' sta]); cd([nwk '_' sta]); % Go to station specific folder to keep things clean . TODO just to cd once. 
