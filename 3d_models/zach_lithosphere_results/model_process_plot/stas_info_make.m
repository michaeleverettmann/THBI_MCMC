%% print some statistics for the LAB and MLD values within the regions
clear all
load('compiled_results_standard.mat');
load('stas_collated.mat');

stas_info = struct('nwk',{mdls.nwk},'sta',{mdls.sta},'lat',mdls.lat,'lon',mdls.lon,'region',which_region,'region_names',regions_str);


% make region 4: NE 
appf = appalachian_front;
gref = grenville_front;

appf_use = appf(appf(:,1)>=-74,:);

for is = 1:length(stas_info.sta)
    if stas_info.region(is) ~=0, continue; end
    if stas_info.lat(is) < 41, continue; end
    if stas_info.lat(is) > 46, continue; end
    if stas_info.lon(is) < -73.76, continue; end    
    if ~geog_inbounds(stas_info.lat(is),stas_info.lon(is)), continue, end
    if stas_info.lon(is) > 0.6+interp1(appf_use(:,2),appf_use(:,1),stas_info.lat(is),'nearest','extrap') 
        stas_info.region(is) = 4;
    end
end
stas_info.region_names(4)="nea";


figure(90),clf,hold on
geogax(gca);
plot(stas_info.lon,stas_info.lat,'^k');
cols = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5];
for ir = 1:4
    plot(stas_info.lon(stas_info.region==ir),stas_info.lat(stas_info.region==ir),'^','markerfacecolor',cols(ir,:))
end
plot(appf(:,1),appf(:,2),'m','linewidth',2)

save('stas_info.mat','stas_info');
