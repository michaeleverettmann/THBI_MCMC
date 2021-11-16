%% Initialise
% proj = struct('name','US');
% proj.dir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/US'; %bb2021.09.27 - This info is going in run_all_stas.m now

error('bb2021.11.15: the proj.mat file is currently generated in evdata1_database.m. Run that, not project_datails.m, or else you will miss important paths. Can also run RUN_prep_data.m.')

%% directories 
proj.rawdatadir = [paths.THBIpath '/data/STAsrawdat/' ]; '/IDKYETVolumes/data/THBI/US/STAsrawdat'; % bb2021.09.27 /volumes/data exists after mounting the NAS in the department.
proj.STAinversions = [paths.THBIpath '/data/STAsinv'];

%% Add matguts to load data to the path
addpath([proj.dir,'/matguts']);
