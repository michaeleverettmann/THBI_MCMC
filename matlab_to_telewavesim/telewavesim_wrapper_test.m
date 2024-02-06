clear; 
pyenv('Version', '~/opt/anaconda3/envs/tws/bin/python', ... % Use anaconda environment where telewavesim is installed. This affects the entire Matlab session. 
    'ExecutionMode','OutOfProcess'); % ERROR ALERT Could not import numpy if using an anconda environment. Matlab would simply crash. However, setting executionMode=OutOfProcess fixed that for me. https://www.mathworks.com/matlabcentral/answers/502458-why-can-py-numpy-array-not-be-resolved

% Simple test parameters. 
modfile = './demo.txt'; % Was created by writeTELEWAVE...m
wvtype = 'P';
npts = 3000;
dt = 0.01;
dp = 2000.0;
use_obs = true;
c = 1.5;
rhof = 1027.0;
slow = 0.06;
baz = 0.0;

% Call the Python function
output = py.run_telewavesim_py.run_telewavesim( ...
    modfile, wvtype, npts, dt, dp, use_obs, c, rhof, slow, baz);
output = cell(output); 
traces = double(output{1})'; 
tt = double(output{2})'; 

figure(1); clf; hold on; 
plot(tt, traces); 