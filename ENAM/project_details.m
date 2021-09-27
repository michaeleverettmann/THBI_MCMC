%% Initialise
% proj = struct('name','US');
% proj.dir = '/Users/zeilon/Documents/MATLAB/BayesianJointInv/US'; %bb2021.09.27 - This info is going in run_all_stas.m now

%% directories 
proj.rawdatadir = '/IDKYETVolumes/data/THBI/US/STAsrawdat'; % bb2021.09.27 /volumes/data exists after mounting the NAS in the department.
proj.STAinversions = '/IDKYETVolumes/data/THBI/US/STAsinv';

%% Add matguts to load data to the path
addpath([proj.dir,'/matguts']);
