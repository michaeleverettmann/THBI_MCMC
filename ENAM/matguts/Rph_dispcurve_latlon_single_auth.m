function [ periods, phV_period, SW_Ray_phV ] = Rph_dispcurve_latlon_single_auth(ilat,ilon,options)
    arguments
        ilat
        ilon
        options.dataset = []
    end
    
%  Function to obtain a composite dispersion curve at any lat/lon point by
%  combination of Earthquake and Ambient noise dispersion curves. Each
%  dispersion curve is interpolated from individual phase velocity maps at
%  the frequencies of interest. If the requested lat/lon point is outside
%  the grid of data, this function will return phV_freq = nan.

paths = getPaths(); 
seismoddir = [paths.models_seismic '/']; 

%% EQ data:
if strcmp(options.dataset, 'DALTON') ; % Dalton et al., 2019? 
    ddir = [seismoddir,'US_RAYLEIGH_EQ_phV_DALTON/'];
    periods = [25,40,50,60,80,100,120,140,180]';
    phV_period = disp_curve_EQ_latlon(periods,ilat,ilon,ddir);
    [periods,iE] = sort(periods);
    phV_period = phV_period(iE);
elseif strcmp(options.dataset, 'SHEN_ANT'); % Shen and Ritzwoller, 2016
    datadir = [seismoddir,'US_RAYLEIGH_ANT_phV_SHEN/'];
    ANperiods = [8:2:32,36,40]';
    ANphV_period = disp_curve_AN_latlon(ANperiods,ilat,ilon,datadir);
elseif strcmp(options.dataset, 'EKSTROM'); 
    % AN data - instead, try Ekstrom. 
    datadir = [seismoddir,'US_RAYLEIGH_ANT_phV_EKSTROM/'];
    periods = [5, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40]';
    phV_period = disp_curve_AN_latlon(periods,ilat,ilon,datadir);
elseif strcmp(options.dataset, 'LYNNER_EQ_EIK') % Colton's Eikonal earthquake phase velocities. 
    warning('No colton eq eikonal yet'); 
elseif strcmp(options.dataset, 'LYNNER_EQ_HELM') % Colton's Eikonal earthquake phase velocities. 
    warning('No colton EQ_HELM'); 
elseif strcmp(options.dataset, 'LYNNER_ANT') % Colton's Eikonal earthquake phase velocities. 
    warning('No colton ANT'); 
elseif strcmp(options.dataset, 'JOSH_ANT') % Colton's Eikonal earthquake phase velocities. 
    warning('No josh ANT'); 
elseif strcmp(options.dataset, 'composite'); 
    warning('Bring in the "composite" stuff from original Rph_dispcurve_latlon.m'); 
end

%% delete any nans
kill = isnan(phV_period);
periods(kill) = [];
phV_period(kill) = [];


if ~isempty(phV_period)
    [periods,iT] = sort(periods);
    phV_period = phV_period(iT);
    SW_Ray_phV = struct('periods',periods,'phV',phV_period,'sigma',[]);
else
    SW_Ray_phV=[];
end

%% Plot. Keep this commented usually
% figure(189); clf; hold on; 
% plot(phVs, periods)
% xlabel('Phase velocity composite')
% ylabel('Period')

end

