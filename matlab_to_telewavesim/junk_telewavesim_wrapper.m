function [traces, tt] = telewavesim_wrapper(modfile, ...
    slow, wvtype, dt, use_obs, dp); 
clear; 
pyenv('Version', '~/opt/anaconda3/envs/tws/bin/python', ... % Use anaconda environment where telewavesim is installed. This affects the entire Matlab session. 
    'ExecutionMode','OutOfProcess'); % ERROR ALERT Could not import numpy if using an anconda environment. Matlab would simply crash. However, setting executionMode=OutOfProcess fixed that for me. https://www.mathworks.com/matlabcentral/answers/502458-why-can-py-numpy-array-not-be-resolved

% Simple test parameters. 
modfile = './demo.txt'; % Was created by writeTELEWAVE...m
wvtype = wvtype(1); % wvtype = 'P';
npts = 3000;
% dt = 0.01; % Sampling period
% dp = 2000.0;
% use_obs = true; % Whether to use OBS or not
c = 1.5; % Water speed
rhof = 1027.0; % Water density
% slow = 0.06; % Ray parameter % slow = 0.06; 
baz = 0.0;

% Call the Python function
output = py.run_telewavesim_py.run_telewavesim( ...
    modfile, wvtype, npts, dt, dp, use_obs, c, rhof, slow, baz);
output = cell(output); 
traces = double(output{1})'; 
tt = double(output{2})'; 

figure(1); clf; hold on; 
plot(tt, traces); 
end