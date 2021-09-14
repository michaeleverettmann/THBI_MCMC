%% script to put various things on path neded for inversion 

% turn off warning about name conflict with matlab builtin isstring
warning('off','MATLAB:dispatcher:nameConflict');
% set time zone:
setenv('TZ','America/Los_Angeles')

basedir = '/Users/brennanbrunsvik/MATLAB/'; % BB Confident in this path % '/Users/zeilon/Documents/MATLAB/'; % change these to your values
bayesdir = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/'; % '/Users/zeilon/Documents/MATLAB/BayesianJointInv/'; % 2021.08.04_TODO % change these to your values

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
addpath('/Users/brennanbrunsvik/MATLAB/fastBSpline'); %  '/Users/zeilon/Dropbox/MATLAB/lib/fastBSpline'); % google "fastBSspline MATLAB" for this
% path to propagator matrix running dir.
addpath([bayesdir,'matlab_to_propmat']); % provided
% path to mineos running dir.
addpath([bayesdir,'matlab_to_mineos']); % bb2021.08.04 changed basedir to bayesdir. Pretty sure basedir was a typo. % provided
% path to seiz models
addpath('/Users/brennanbrunsvik/MATLAB/seizmo/models'); %'~/Dropbox/MATLAB/lib/seizmo/models') % Google "MATLAB seizmo" for this
% path to HV kernels functions
% addpath('~/Work/codes/HV_Tanimoto/matlab_to_HVkernel') % BB commented out 2021.08.04 % only if doing ellipticity ratios
% path to Rayleigh wave dispersion curve dir
addpath('/Users/brennanbrunsvik/MATLAB/seis_tools-master/surface_waves'); % '/Users/zeilon/Dropbox/MATLAB/seis_tools/surface_waves'); % on my github
% path to gaussfit dir
% addpath('~/Dropbox/MATLAB/lib/gaussfit'); % BB commented 2021.08.04 % don't think used any more

addpath('/Users/brennanbrunsvik/MATLAB/seis_tools-master'); % bb2021.08.04 Adding because some functions were missing that are here. I'm worried this might break something though. Why was it not already added? 
% addpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup_alt'); % bb2021.08.04 Need taup

addpath('/Users/brennanbrunsvik/MATLAB/EilonmyFUNCTIONS'); % bb2021.08.04 added to get "maxab"

% Some things needed for TauP
javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar'); 
javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/TauP-2.1.1.jar');
addpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup'); % This might break shit. It's Seizmo. 
addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/misc/for_seizmo'); % bb2021.08.04 Just taking the absolutely necessary Seizmo tools. 


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