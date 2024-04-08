run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

mdls = load(fresults).mdls; % Get station positions and uncertainties. 

lat = mdls.lat; 
lon = mdls.lon; 

param_plot = ["VSsig1"]; 
i_param = 1; 
param_i = param_plot(i_param); 


n_stations = size(mdls.lat,1); 
parm_all = nan(n_stations, 1); 

% Gather data from each station
for i_station = 1:n_stations; 
    model =  mdls.model{i_station}; 
    parm = model.(param_i); 
    
    idepth = 20; 
    parm = parm(idepth,:); 
    parm = parm(2) - parm(1); 
    parm_all(i_station) = parm;
    
end; 

%%
figure(1); clf; hold on; 
scatter(lon, lat, 20, parm_all, 'filled'); 
colorbar(); 