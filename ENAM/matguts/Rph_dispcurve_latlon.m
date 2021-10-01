function [ periods, phVs ] = Rph_dispcurve_latlon(ilat,ilon,transT)
% [ periods, phVs ] = Rph_dispcurve_latlon(ilat,ilon,transT)
% 
%  Function to obtain a composite dispersion curve at any lat/lon point by
%  combination of Earthquake and Ambient noise dispersion curves. Each
%  dispersion curve is interpolated from individual phase velocity maps at
%  the frequencies of interest. If the requested lat/lon point is outside
%  the grid of data, this function will return phV_freq = nan.

EQauth = 'DALTON'; 
ANauth = 'EKSTROM'; 

if nargin<3 || isempty(transT)
    transT = 33;
end

paths = getPaths(); 
seismoddir = [paths.models_seismic '/']; 

% 
% seismoddir = '/Volumes/data/models_seismic/';
% if ~exist(seismoddir,'dir')
%     try
%         seismoddir = '/Volumes/eilon_data/models_seismic/';
%     catch
%         error('NO SEISMOD DIR FOUND');
%     end
% end

%% EQ data:
if strcmp(EQauth, 'DALTON'); 
    ddir = [seismoddir,'US_RAYLEIGH_EQ_phV_DALTON/'];
    Eperiods = [25,40,50,60,80,100,120,140,180]';
    EphV_period = disp_curve_EQ_latlon(Eperiods,ilat,ilon,ddir);
    [Eperiods,iE] = sort(Eperiods);
    EphV_period = EphV_period(iE);
else
    error('bb2021.09.30 I didnt build code for this earthquake phase velocity author yet. Should be easy to do this.') 
end


%% AN data:
if strcmp(ANauth, 'SHEN'); 
    % Shen and Ritzwoller, 2016
    datadir = [seismoddir,'US_RAYLEIGH_ANT_phV_SHEN/'];
    ANperiods = [8:2:32,36,40]';
    ANphV_period = disp_curve_AN_latlon(ANperiods,ilat,ilon,datadir);
elseif strcmp(ANauth, 'EKSTROM'); 
    % AN data - instead, try Ekstrom. 
    datadir = [seismoddir,'US_RAYLEIGH_ANT_phV_EKSTROM/'];
    ANperiods = [5, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40]';
    ANphV_period = disp_curve_AN_latlon(ANperiods,ilat,ilon,datadir);
else
    error('bb2021.09.30 I didnt build code for this ambient noise phase velocity author yet. Should be easy to do this.')
end 


%% composite
p1 = ANperiods(ANperiods<min(Eperiods));
v1 = ANphV_period(ANperiods<min(Eperiods));

p3 = Eperiods(Eperiods>transT);
v3 = EphV_period(Eperiods>transT);

p2 = [min(Eperiods):2:transT]';
v2 = mean([interp1(Eperiods,EphV_period,p2),interp1(ANperiods,ANphV_period,p2)],2);

periods = [p1;p2;p3];
phVs    = [v1;v2;v3];

%% delete any nans
kill = isnan(phVs);
periods(kill) = [];
phVs(kill) = [];

%% Plot. Keep this commented usually
% figure(189); clf; hold on; 
% plot(phVs, periods)
% xlabel('Phase velocity composite')
% ylabel('Period')

end

