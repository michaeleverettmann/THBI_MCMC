%% script to put various things on path neded for inversion 
% Because this is so project dependent, probably will need to put in in the proj.name folder eventually. 
proj = struct('name', 'ENAM'); % Subfolder where this inversions stuff is stored. 

% turn off warning about name conflict with matlab builtin isstring
warning('off','MATLAB:dispatcher:nameConflict');
% set time zone:
setenv('TZ','America/Los_Angeles')

% Get home directory. Everything should be placed about this. 
% Try to keep paths relative to home the same on EVERY computer. 
[~,hd] = system('echo ~'); % TODO use getenv('HOME') to get home directory
hd = strip(hd); 

basedir = [hd '/MATLAB/']; % BB Confident in this path % '/Users/zeilon/Documents/MATLAB/'; % change these to your values
bayesdir = [hd '/Documents/UCSB/ENAM/THBI_ENAM/']; % '/Users/zeilon/Documents/MATLAB/BayesianJointInv/'; % 2021.08.04_TODO % change these to your values

% You'll likely need to get all of these functions - email me if not
% obvious. Some turn out not to be used any more, so comment out and then
% let me know if they come up and you need them. 

% path of all inversion main functions
addpath(bayesdir)
% path of all inversion sub functions
addpath([bayesdir,'functions']) % inside this folder
% path of sub functions that I got online but don't know if I can legally
% distribute. bb2021.09.14
addpath(genpath([bayesdir, 'functions/functionsExternal'])); 
% path to fast spline func.
% Should go to this folder and run the CompileMex something file. 
addpath([hd '/MATLAB/fastBSpline']); %  '/Users/zeilon/Dropbox/MATLAB/lib/fastBSpline'); % google "fastBSspline MATLAB" for this
% path to propagator matrix running dir.
addpath([bayesdir,'matlab_to_propmat']); % provided
% path to mineos running dir.
addpath([bayesdir,'matlab_to_mineos']); % bb2021.08.04 changed basedir to bayesdir. Pretty sure basedir was a typo. % provided
addpath([bayesdir,'matlab_to_hv_kernel']); 
% path to seiz models
addpath([hd '/MATLAB/seizmo/models']); %'~/Dropbox/MATLAB/lib/seizmo/models') % Google "MATLAB seizmo" for this
% path to HV kernels functions
% addpath('~/Work/codes/HV_Tanimoto/matlab_to_HVkernel') % BB commented out 2021.08.04 % only if doing ellipticity ratios
% path to Rayleigh wave dispersion curve dir
addpath([hd '/MATLAB/seis_tools-master/surface_waves']); % '/Users/zeilon/Dropbox/MATLAB/seis_tools/surface_waves'); % on my github
% path to gaussfit dir
% addpath('~/Dropbox/MATLAB/lib/gaussfit'); % BB commented 2021.08.04 % don't think used any more

addpath([hd '/MATLAB/seis_tools-master']); % bb2021.08.04 Adding because some functions were missing that are here. I'm worried this might break something though. Why was it not already added? 
% addpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup_alt'); % bb2021.08.04 Need taup

addpath([hd '/MATLAB/EilonmyFUNCTIONS']); % bb2021.08.04 added to get "maxab"

% Some things needed for TauP
javaaddpath([hd '/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar']); 
javaaddpath([hd '/MATLAB/seizmo/mattaup/lib/TauP-2.1.1.jar']);
javaaddpath([hd '/MATLAB/iris/IRIS-WS-2.0.18.jar']); % bb2021.09.27 moving this path from evdata1_database.m. Should be just for downloading station data from IRIS. 

addpath([hd '/MATLAB/seizmo/mattaup']); % This might break things. It's Seizmo. 
addpath([hd '/Documents/UCSB/ENAM/THBI_ENAM/functions/misc/for_seizmo']); % bb2021.08.04 Just taking the absolutely necessary Seizmo tools. 


% addpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup'); % bb2021.08.04 Need taup. 
% addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/misc/for_seizmo'); % bb2021.08.04 Just taking the absolutely necessary Seizmo tools. 
% % javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar'); % Needed for taup I think.
% javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar'); % Needed for taup I think.
% javaaddpath('/Users/brennanbrunsvik/MATLAB/TauP-master/src/main/java/edu'); % Needed for taup I think
% javaaddpath('/Users/brennanbrunsvik/MATLAB/TauP-master/src/main/java/edu/sc/seis/TauP/SphericalCoords.java')
% javaaddpath('/Users/brennanbrunsvik/MATLAB/TauP-master/gradle/wrapper/gradle-wrapper.jar')
% javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/TauP-2.1.1.jar'); 

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict');

global prem_anisotropic prem_isotropic
prem_anisotropic = prem_perfect('SPVW',0.25);
prem_isotropic = prem;

% Some special paths which you list here once. Using this set and get paths
% approach means we don't have to use global variables, which work poorly
% (if at all) with parfor. 
setPaths(hd,proj); 
% paths = getPaths(); 
function setPaths(hd,proj); 
% system('FILE=/usr/local/bin/timeout;if test -f "$FILE"; then; echo "$FILE exists.";fi')
if isfile('/usr/local/bin/timeout'); % Would be better to do system('which timeout'), but system does not have the normal $PATH, and cannot always find timeout. 
    timeoutpath = '/usr/local/bin/timeout'; 
elseif isfile('/usr/bin/timeout'); 
    timeoutpath = '/usr/bin/timeout'; 
else
    error('Brennan warning: cannot find timeout path on your computer. Needed for Mineos. Find it in shell script using which timeout or which gtimeout, then enter it as timeoutpath.')
end

paths = struct('CADMINEOS', [hd '/Documents/repositories/Peoples_codes/CADMINEOS'],...
                   'PropMat', [hd '/Documents/repositories/Peoples_codes/PropMat'],...
                   'HV_ellipticity', [hd '/Documents/repositories/Peoples_codes/HV_ellipticity'],...
                   'timeout', timeoutpath,...
                   'THBIpath', [hd '/Documents/UCSB/ENAM/THBI_ENAM'],...
                   'rawdatadir', ['/Volumes/data/',proj.name,'/THBI/STAsrawdat/'],...
                   'STAinversions', ['/Volumes/data/',proj.name,'/THBI/STASinv/'],...
                   'models_seismic', [hd '/Documents/repositories/data/models_seismic']); 
paths.rawdatadir = [paths.THBIpath '/data/STAsrawdat/']; % bb2021.09.28 if this takes too much local storage, put it on external drives. 
paths.STAinversions = [paths.THBIpath '/data/STASinv/']; % bb2021.09.28 if this takes too much local storage, put it on external drives. 
% 'rawdatadir', ['/Volumes/data/',proj.name,'/THBI/STAsrawdat/'],... % bb2021.09.28 Could use these locations to keep data on NAS. Not so important though if you aren't downloading a bunch of waveform data I think. 
% 'STAinversions', ['/Volumes/data/',proj.name,'/THBI/STASinv/'] ); 
% save([paths.THBIpath '/misc/paths.mat']); 
matlab.io.saveVariablesToScript('pathsAutoGen.m', 'paths')

addpath(paths.models_seismic); % bb2021.11.02 Needed for some functions which access data. I thought this was already added somewhere but I don't see where. 
end
