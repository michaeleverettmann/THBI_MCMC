function [] = fun_mapbackground(lon, lat, dat, this_title, lolim, lalim, options); 
    arguments
        lon
        lat
        dat
        this_title
        lolim = []
        lalim = []
        options.stanames = []
    end
% Make the maps background for a few other scripts. 

lolim = [-87, -71]; 
lalim = [27 ,  46]; 

nexttile; cla; hold on;
m_proj('miller', 'long',lolim + [-2 2],'lat',lalim + [-2 2]);
m_coast('patch',[1 1 1]); 

[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);

[stax, stay] = m_ll2xy(lon, lat);  
scatter(stax, stay, 40, dat, 'filled'); 
title(this_title, 'Fontweight', 'normal'); 

colorbar(); 

turbo_arr = turbo(); 
colormap(turbo_arr(end:-1:1,:)); 

if ~isempty(options.stanames); 
    for ista = 1:length(lon); 
        text(stax, stay, options.stanames); 
    end
end

end