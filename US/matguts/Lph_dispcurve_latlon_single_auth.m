function [ periods, phV_period, SW_Lov_phV ...
    ] = Lph_dispcurve_latlon_single_auth(ilat,ilon,options)
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

%% Love data
if strcmp(options.dataset, 'SW_Lov_phV'); 
    error('We now have multiple Love wave phase velocity datasets. Specify which one. brb20240703')
end 

if strcmp(options.dataset, 'SW_Lov_phV_har') % Colton's Eikonal earthquake phase velocities. 
    ddir = [seismoddir,'US_LOVE_HARIHARAM/'];
    periods = [35, 40, 45, 50, 60, 75]';
    [periods, phV_period] = load_love_hariharam(ilat, ilon, ddir, 'periods', periods);
elseif strcmp(options.dataset, 'SW_Lov_phV_eks'); 
    [periods,phV_period]  = Lph_dispcurve_latlon( ilat,ilon ); % Ekstrom dataset. 
end

%% delete any nans
kill = isnan(phV_period);
periods(kill) = [];
phV_period(kill) = [];


if ~isempty(phV_period)
    [periods,iT] = sort(periods);
    phV_period = phV_period(iT);
    SW_Lov_phV = struct('periods',periods,'phV',phV_period,'sigma',[]);
else
    SW_Lov_phV=[];
end

%% Plot. Keep this commented usually
% figure(189); clf; hold on; 
% plot(phVs, periods)
% xlabel('Phase velocity composite')
% ylabel('Period')

end

