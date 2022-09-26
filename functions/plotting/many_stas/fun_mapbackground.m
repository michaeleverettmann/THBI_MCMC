function [] = fun_mapbackground(lon, lat, dat, this_title, lolim, lalim); 

lolim = [-87, -74]; 
lalim = [30 ,  43]; 

nexttile; cla; hold on;
m_proj('lambert', 'long',lolim + [-2 2],'lat',lalim + [-2 2]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);

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

end