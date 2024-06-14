function [maxHK] = hk_maximum_possible_value(RF,tt); 
% RF: receiver functions. number of time steps by number of receiver
% functions. 
% tt: time. 

% Super simplified way to estimate the maximum value that an HK stack could
% ever achieve for some receiver functions. 

% predata.HKstack_P.waves.rf
% predata.HKstack_P.waves.tt

tmin = 2; % Minimum time to consider in RFs. Width of pulse? 
tmax = 40; % Max time 

maxHK = 0; 
for iWave = [1:size(RF,2)]; 

    inWind = and(tt>=tmin,tt<=tmax); 
    [pks,locs] = findpeaks( RF(inWind),tt(inWind) ); 
    [pksSort,indSort] = sort(pks); 
    pksSort = flip(pksSort); 

    % Maximum obtainable HK value. 
    % Assume Ps phase is the largest, PpPs is second largest, PpSs,
    % PsPs are the most negative value. 
    maxHKStack = 0.5 * pksSort(1) + 0.3 * pksSort(2) - 0.2 * pksSort(end); 
    maxHK = maxHK + maxHKStack; % Simply sum across all HK stacks. 

end