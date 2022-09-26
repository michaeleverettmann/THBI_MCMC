%% Script to establish a database of stations and events for body wave study
% Don't use clear at the beginning of this script. Sometimes I run it after
% already having provided network/station name in RUN_one_sta.m
ifrunall = true;
ifdownloadseed = false; 

run('../a0_STARTUP_BAYES.m')
paths = getPaths(); 
proj = struct('name','ENAM');
proj.dir = [paths.THBIpath '/' proj.name];

%% Station parameters
sta_latlims = [21 48]; % [min_lat max_lat] for stations
sta_lonlims = [-93, -63.4]; % [min_lon max_lon] for stations
sta_chans = 'BH*,HH*'; % channel codes to search for
starttime = '1970-01-01 00:00:00';
startbytime = '2030-01-01 00:00:00';
min_longevity_yrs = 0.6;

%% Event parameters
mag_lims = [5.0 7.8];
dep_lims = [0 1000]; % set to [0 1000] by default (km)
gc_lims  = [30 75];
% startafter = '1990-01-01 00:00:00'; % earliest evtime (yyyy-mm-dd HH:MM:SS)

%% Data parameters
phases = {'P','S'};
samprate = 40;
datwind =  [-100 100];

%% ID for IRIS DMC request
IRIS_ID = 'bbrunsvik';


%% GET TO WORK
wd = pwd;
addpath('matguts');

%% Make directory structure
% main database directory
proj.dir = regexprep(proj.dir,'~',getenv('HOME'));
if ~strcmp(proj.dir(end),'/'),proj.dir = [proj.dir,'/']; end
if exist(proj.dir,'dir')~=7, mkdir(proj.dir); end
% info files directory 
proj.infodir = [proj.dir,'INFO/'];
if exist(proj.infodir,'dir')~=7, mkdir(proj.infodir); end
% response files directory 
proj.respdir = [proj.dir,'INFO/RESP/'];
if exist(proj.respdir,'dir')~=7, mkdir(proj.respdir); end
% data files directory 
proj.rawdatadir = paths.rawdatadir;
proj.STAinversions = paths.STAinversions;
if exist(proj.rawdatadir,'dir')~=7, mkdir(proj.rawdatadir); end
if exist(proj.STAinversions,'dir')~=7, mkdir(proj.STAinversions); end

% save project details
save([proj.dir,'project_details.mat'],'proj');

%% Add matguts to load data to the path
addpath([proj.dir,'/matguts']);


cd(proj.dir)

if ~ifrunall
    return 
end

%% Write request details information
request_details_all = struct(...
        'phases',	{phases},... % need double {{ }}
        'gclims',   gc_lims,...
        'maglims',  mag_lims,...
        'samprate', samprate,...
        'datwind',  datwind);

save([proj.infodir,'/data_request_details.mat'],'request_details_all');

if ~ifrunall
    return
end

%% Get station + channel information
% bb2021.09.27 Make sure you have IRIS-WS-2.0.18.jar or something similar added to java path for this part. Should have been added to the path in a0_STARTUP_BAYES.M
% save station request info
stations_request = struct('lat_lims',sta_latlims,'lon_lims',sta_lonlims,'chans',sta_chans,...
                          'starttime',starttime,'startbytime',startbytime,'min_longevity_yrs',min_longevity_yrs);

% grab stations
stations_IRIS = irisFetch.Stations('station','*','*','*','BH?',...
    'boxcoordinates',[sta_latlims,sta_lonlims],'StartAfter',starttime,'StartBefore',startbytime);

%%% For loading station details later. Sort of a hack...
for ista = 1:length({stations_IRIS.NetworkCode}); 
    stadeets_netsta = stations_IRIS(ista); 
    save(sprintf('%sstadeets_%s_%s.mat',...
        paths.rawdatadir,stadeets_netsta.NetworkCode,stadeets_netsta.StationCode),...
        'stadeets_netsta'); 
end
%%%

%only include stations satisfying longevity
starter = datenum({stations_IRIS.StartDate}');
for ed = 1:length(stations_IRIS); % This isn't vectorized because it's a cell function and that was a pain. bb2021.09.27
    if length(stations_IRIS(ed).EndDate) < 3; % We get [] for active stations. I chose length < 3 arbitrarily. Even 1 should work, but 3 feels safer.  
        stations_IRIS(ed).EndDate = '2099-12-31 00:00:00.000'; 
    end
end
ender = datenum({stations_IRIS.EndDate}'); ender(ender>now)= now;
stations_IRIS = stations_IRIS((ender - starter)/365.25 >= min_longevity_yrs)


stainfo = struct('stas',{{stations_IRIS.StationCode}'},...
                 'nwk',{{stations_IRIS.NetworkCode}'},...
                 'slats',[stations_IRIS.Latitude]',...
                 'slons',[stations_IRIS.Longitude]',...
                 'selevs',[stations_IRIS.Elevation]',...
                 'ondate',datenum({stations_IRIS.StartDate}'),...
                 'offdate',datenum({stations_IRIS.EndDate}'),...
                 'ondate_str',{{stations_IRIS.StartDate}'},...
                 'offdate_str',{{stations_IRIS.EndDate}'},...
                 'nstas',length(stations_IRIS));  
             
[stainfo] = stainfo_unique(stainfo);

% parse channels             
chans = cell(stainfo.nstas,3);
chandips = nan(stainfo.nstas,3);
chanazs = nan(stainfo.nstas,3);
nchans = zeros(stainfo.nstas,1);
for is = 1:stainfo.nstas
    nchan = length(stations_IRIS(is).Channels);
    tempchans = cell(1,nchan);
    tempdips = nan(1,nchan);
    tempazs = nan(1,nchan);
    for ic = 1:nchan
        tempchans(ic) = {stations_IRIS(is).Channels(ic).ChannelCode};
        tempdips(ic) = stations_IRIS(is).Channels(ic).Dip;
        tempazs(ic) = stations_IRIS(is).Channels(ic).Azimuth;
    end
    [stachans,indch] = unique(tempchans);
    nchans(is) = length(stachans);
    chans(is,1:nchans(is)) = stachans;
    chandips(is,1:nchans(is)) = tempdips(indch);
    chanazs(is,1:nchans(is)) = tempazs(indch);
end
chandips(cellfun('isempty',chans)) = nan;
chanazs(cellfun('isempty',chans)) = nan;
stainfo.nchans = nchans;
stainfo.chans = chans;
stainfo.chandips = chandips;
stainfo.chanazs = chanazs;

% % % brb2022.02.28 Get receiver function traces downloading. 
% % paths = getPaths(); 
% % thisDir = pwd(); 
% % rfWaveformFolder = [paths.models_seismic, '/US_EARS/IRIS_EARS']; 
% % cd(rfWaveformFolder)
% % [inf1, inf2] = system(['conda activate seis_gen';...
% %     'python getEarsRFByStation.py';...
% %     'python getEarsRFByStation_unpack.py']); % Run python codes to download IRIS data. 
% % cd(thisDir); 
% % % End getting receiver function waveforms. 

save([proj.infodir,'stations'],'stainfo','stations_IRIS','stations_request');
save([proj.infodir 'proj.mat'], 'proj'); 

if ~ifrunall
    return
end

if ifdownloadseed; 
    %% Response SAC_PZ files
    % Build and send BREQFAST request file for dataless seed                     
    addpath('~/Dropbox/MATLAB/seis_tools/breqfasting/');
    breq_fast_request([proj.name,'_dataless'],IRIS_ID,{stations_IRIS.StationCode}','*',{stations_IRIS.NetworkCode}','',{stations_IRIS.StartDate}',{stations_IRIS.EndDate}','dataless_SEED',[proj.name,'_dataless_request']);
    movefile([proj.name,'_dataless_request'],proj.respdir);
    return 
    %% ======= WAIT FOR NOTIFICATION THAT DATALESS SEED IS ON SERVER ======= %%
    % Download and process dataless seed  
    breq_fast_dataless_PZprocess([proj.name,'_dataless'],IRIS_ID,proj.respdir,{'BH*','HH*'},1 )     
    return
end
                        