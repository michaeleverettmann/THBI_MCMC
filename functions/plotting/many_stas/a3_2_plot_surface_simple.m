function [] = a3_2_plot_surface_simple(llminmax, options)
    arguments
        llminmax 
        options.stax = []
        options.stay = []
        options.stav = []
        options.stalon = [] 
        options.stalat = [] 
        options.xgrid = []
        options.ygrid = [] 
        options.vgrid = []
        options.fignum = 1; 
        options.title = []; 
    end

% llminmax = [lonmin, lonmax, latmin, latmax]. Is used to define map
% projection and relation to x y coordinates. 

lonmin = llminmax(1); 
lonmax = llminmax(2); 
latmin = llminmax(3); 
latmax = llminmax(4); 

%% Establish projection

figure(options.fignum); clf; hold on; 
set(gcf, 'color', 'white'); 
m_proj('lambert', 'long',[lonmin, lonmax],'lat',[latmin, latmax]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);


%% Plot surface
if ~isempty(options.xgrid); 

    if ~isempty(options.stax); 
        pt_dist_nan = 0.02; 
        pt_dist = zeros(size(options.xgrid)); 
        for ipt = 1:(size(options.xgrid, 1) * size(options.xgrid,2)); 
            pt_dist(ipt) = min(sqrt((options.xgrid(ipt) - options.stax).^2 + ...
                 (options.ygrid(ipt) - options.stay).^2)); 
        end
        options.vgrid(pt_dist > pt_dist_nan) = nan; 
    end

    contourf(options.xgrid, options.ygrid, options.vgrid, 15,...
        'LineStyle','none'); 
%     m_contourf(options.xgrid, options.ygrid, options.vgrid, 15); 
end


[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);
m_grid('box','fancy','linestyle','-','gridcolor',.5 .*[1,1,1],'backcolor',[.3 .75 1]);
%%% TODO find a way to nan all the data that leaks out of the projection
%%% box. 


if ~isempty(options.title); 
    title(options.title, 'fontweight', 'normal')
end





%% Handle choices in terms of plotting stations from lon/lat, x/y, or not plotting stations. 
if isempty(options.stax); 
    if ~isempty(options.stalon); 
        [stax, stay] = m_ll2xy(options.stalon, options.stalat);  
    else;
        stax = []; 
        stay = []; 
    end
else
    stax = options.stax; 
    stay = options.stay; 
end
if isempty(options.stav); 
    options.stav = ones(size(options.stax)); 
end

%% Plot stations. 
if ~isempty(stax);     
    scatter(stax, stay, 30, options.stav, 'filled', 'MarkerEdgeColor','k'); 
end

colorbar(); 
cmap = turbo(); 
cmap = cmap(end:-1:1,:); 
colormap(cmap); 
try
    caxis([min(options.stav), max(options.stav)]);
catch error_rep
    fprintf('%s', getReport(error_rep)); 
end

end