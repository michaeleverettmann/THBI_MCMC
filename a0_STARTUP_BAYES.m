%% This script mostly handles your paths to other repositories or folders. 
% You will likely need to read through the required files and find them. 
% Most paths that need manual changing should be here. 
% Don't use clear at the beginning of this script. Sometimes I run it after already having provided network/station name in RUN_one_sta.m

proj = struct('name', 'ENAM'); % Subfolder that can contain things for your application of MCMC to your own data. TODO not sure this matters now. 

setenv('TZ','America/Los_Angeles'); % set time zone. Maybe not needed? 

% Get home directory. Try to keep path structures, relative to home, the same on any computer. 
[~,hd] = system('echo ~'); % TODO use getenv('HOME') to get home directory. Need to test on linux. 
hd = strip(hd); 

basedir = [hd '/MATLAB/']; % Contains a few matlab general directories. 
bayesdir = [hd '/Documents/UCSB/ENAM/THBI_ENAM/']; % Where your THBI_MCMC folder (github repo) stored

%% MCMC provided subfunctions
% You'll likely need to get all of these functions - email me if not obvious. Some turn out not to be used any more, so comment out and then let me know if they come up and you need them. 
% path of all inversion main functions
addpath(bayesdir)
% path of all inversion sub functions
addpath([bayesdir,'functions']) % inside this folder
addpath([bayesdir,'functions/hKappaFunk']) % For doing some things with h-kappa stacks. 
addpath([bayesdir,'functions/plotting/progress']) % For doing some things with h-kappa stacks. 
addpath([bayesdir,'functions/master_subfuncs']); % Trying to break up master par into many subfunctions
addpath([bayesdir,'functions/handle_errors']);
addpath([bayesdir,'functions/plotting']); 
addpath([bayesdir,'functions/plotting/plot_data_heatmap']); 
addpath([bayesdir,'functions/plotting/sensitivity']); 
addpath([bayesdir,'functions/plotting/receiver_functions']); 


%% Path of sub functions that I got online but don't know if I can distribute. bb2021.09.14
addpath(genpath([bayesdir, 'functions/functionsExternal'])); 
addpath([hd '/MATLAB/fastBSpline']); % https://www.mathworks.com/matlabcentral/fileexchange/32509-fast-b-spline-class. Path to fast spline func. Should go to this folder and run the CompileMex something file. google "fastBSspline MATLAB" for this. 
addpath([bayesdir,'matlab_to_propmat']); % path to propagator matrix running dir. provided

%%% brb20240621 some stuff for changing the mineos version
% addpath([bayesdir,'matlab_to_mineos']); % path to mineos running dir. provided
disp('brb20240621 Temporarily not adding matlab_to_mineos to path. '); 
addpath([hd '/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_mineos_vJBR/functions']); % Remove when ready
addpath([hd '/Documents/repositories/Peoples_codes/MINEOS_synthetics/run_MINEOS']); % Move when ready 
%%%

addpath([bayesdir,'matlab_to_hv_kernel']); % matlab wrapper for HV codes. 
addpath([hd '/MATLAB/seizmo/models']); % path to seiz models. Google "MATLAB seizmo" for this
addpath([hd '/MATLAB/seis_tools-master/surface_waves']); % path to Rayleigh wave dispersion curve dir. '/Users/zeilon/Dropbox/MATLAB/seis_tools/surface_waves'); % on Zach's github
addpath([hd '/MATLAB/seis_tools-master']);  
% addpath('~/Dropbox/MATLAB/lib/gaussfit'); % path to gaussfit dir brb20210804 don't think used any more
% addpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup_alt'); % brb2021.08.04 Need taup
addpath([hd '/MATLAB/EilonmyFUNCTIONS']); % bb2021.08.04 added to get "maxab"

%% telewavesim and it's Python environment. 
% brb20240619 This Python environment should only be needed if you are
% using telewavesim in stead of propmat. 
p_tws = [bayesdir,'matlab_to_telewavesim']; 
addpath(p_tws); 
if count(py.sys.path, p_tws) == 0
    insert(py.sys.path, int32(0), p_tws); % Telewavesim depends on Python, so we also need to modify our Python path. 
end
pyenv('Version', '~/opt/anaconda3/envs/tws/bin/python', ... % Use anaconda environment where telewavesim is installed. This affects the entire Matlab session. TODO define this path in somewhere more obvious. 
    'ExecutionMode','OutOfProcess'); % ERROR ALERT Could not import numpy if using an anconda environment. Matlab would simply crash. However, setting executionMode=OutOfProcess fixed that for me. https://www.mathworks.com/matlabcentral/answers/502458-why-can-py-numpy-array-not-be-resolved

%%
% Some things needed for TauP. Taup used with synthetic receiver functions.
javaaddpath([hd '/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar']); 
javaaddpath([hd '/MATLAB/seizmo/mattaup/lib/TauP-2.1.1.jar']);
javaaddpath([hd '/MATLAB/iris/IRIS-WS-2.0.18.jar']); % bb2021.09.27 moving this path from evdata1_database.m. Should be just for downloading station data from IRIS. So you shouldn't need it.  

addpath([hd '/MATLAB/seizmo/mattaup']); % This might break things. It's Seizmo and can conflict with default variable names.  
addpath([hd '/Documents/UCSB/ENAM/THBI_ENAM/functions/misc/for_seizmo']); % bb2021.08.04 Just taking the absolutely necessary Seizmo tools. 

% Any specific colormaps we want. Careful not to add them all or you might overwrite a default variable. 
addpath([hd '/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri']);
addpath([hd '/Documents/repositories/Base_code/colormaps/ScientificColourMaps7/vik']); 

% HK stacks. From https://github.com/brennanbrunsvik/hk_anis
addpath([hd, '/Documents/repositories/hk_anis/elastic']); % Two paths needed for doing anisotropic HK stacks. 
addpath([hd, '/Documents/repositories/hk_anis/hk_calculation']);

%% Various things that might have been needed for taup, but not in the codes current form. 
% addpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup'); % bb2021.08.04 Need taup. 
% addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/misc/for_seizmo'); % bb2021.08.04 Just taking the absolutely necessary Seizmo tools. 
% % javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar'); % Needed for taup I think.
% javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/MatTauP-2.1.1.jar'); % Needed for taup I think.
% javaaddpath('/Users/brennanbrunsvik/MATLAB/TauP-master/src/main/java/edu'); % Needed for taup I think
% javaaddpath('/Users/brennanbrunsvik/MATLAB/TauP-master/src/main/java/edu/sc/seis/TauP/SphericalCoords.java')
% javaaddpath('/Users/brennanbrunsvik/MATLAB/TauP-master/gradle/wrapper/gradle-wrapper.jar')
% javaaddpath('/Users/brennanbrunsvik/MATLAB/seizmo/mattaup/lib/TauP-2.1.1.jar'); 

%%
global prem_anisotropic prem_isotropic
prem_anisotropic = prem_perfect('SPVW',0.25);
prem_isotropic = prem;


%% Below I set some special paths which will need to be retrieved in functions. 
% Instead of storing them in global variables (which doesn't work with
% parfor), we make a matlab script that will load these paths from a .mat
% file. 

ramDrive = make_ram_drive('ramMb', feature('numcores')*768 ); % Make a folder in ram for faster temporary IO. Tested on a Mac and a few Linux computers. 
setPaths(hd,proj,ramDrive); % Use below function to store the paths. 

function setPaths(hd,proj,ramDrive); 

    % brb20240628 For development on a local computer, might want to store
    % "data" on an external drive. If that external drive is found, assume
    % that data is stored here: 
    data_drive = '/Volumes/extDrive/offload/Users/brennanbrunsvik'; 
    if ~ (exist(data_drive)==7); % Data drive not present. 
        data_drive = hd; % Assume data is stored at hd. File structure from point of data_drive and hd should be the same on your local machine, if you store data on an external drive. 
    end 

    % Need shell function timeout. Might be in different places depending on the computer. 
    if isfile('/usr/local/bin/timeout'); % Would be better to do system('which timeout'), but system does not have the normal $PATH, and cannot always find timeout. 
        timeoutpath = '/usr/local/bin/timeout'; 
    elseif isfile('/usr/bin/timeout'); 
        timeoutpath = '/usr/bin/timeout'; 
    else
        error('Brennan warning: cannot find timeout path on your computer. Needed for Mineos. Find it in shell script using which timeout or which gtimeout, then enter it as timeoutpath.')
    end
    
    paths = struct(    'CADMINEOS',      [hd '/Documents/repositories/Peoples_codes/CADMINEOS'],...
                       'PropMat',        [hd '/Documents/repositories/Peoples_codes/PropMat'],...
                       'HV_ellipticity', [hd '/Documents/repositories/Peoples_codes/HV_ellipticity'],...
                       'timeout',         timeoutpath,...
                       'THBIpath',       [hd '/Documents/UCSB/ENAM/THBI_ENAM'],...
                       'execPath',       ramDrive,... % Working directory during inversion. Is on ram, for fast IO. 
                       'models_seismic', [data_drive '/Documents/repositories/data/models_seismic'], ... % Where is your seismic data stored. 
                       'ramDrive',        ramDrive); 
    
    addpath(sprintf('%s/SEMum2_avg_VS',paths.models_seismic)); % brb2022.08.16 Average global velocity models. For synthetic tests.                 
                   
    paths.rawdatadir =    [paths.THBIpath '/data/STAsrawdat/']; % bb2021.09.28 if this takes too much local storage, put it on external drives. 
    paths.STAinversions = [paths.THBIpath '/data/STASinv/'   ]; % bb2021.09.28 if this takes too much local storage, put it on external drives. 
    
    % Save paths. Sometimes running in parallel caused a glitch here, so try a few times if needed (maybe not be needed now). 
    itries = 0; 
    while itries < 10; 
        try 
            pause_time = rand(1)*.1; % Pause to reduce risk of parallel writing to same file. maybe not needed now. 
            % fprintf('Pausing for %1.3f seconds before loading paths file\n',pause_time); 
            pause(pause_time); 
            matlab.io.saveVariablesToScript([paths.THBIpath '/pathsAutoGen.m'], 'paths'); % Save variables to a script which we can call from any function
            paths = getPaths(); 
            itries = inf; 
        catch e
            fprintf('Error reading or saving paths. Will try again 10 times. Error was: \n%s\n',getReport(e)); 
            itries = itries + 1; 
            if itries > 10; 
                error('Failed at paths 10 times. Giving up. '); 
            end
        end
    end
    
    addpath(paths.models_seismic); % bb2021.11.02 Needed for some functions which access data.  
end
