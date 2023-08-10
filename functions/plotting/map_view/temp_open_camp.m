%%% brb2023/08/09 This script determines how Zach's CAMP digitized file
%%% works. 

load('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/CAMP_digitize/CAMP_fromGao2000f1_nan.mat'); 

figure(1); 
clf; hold on; 
scatter(a(:,2), a(:,3), 'filled'); 

brk = find(isnan(a(:,2))); 
brk = [1; brk; size(a,1)]; 

for ibrk = 1:(length(brk)-1); 
    lon = a(brk(ibrk):brk(ibrk+1),2); 
    lat = a(brk(ibrk):brk(ibrk+1),3); 
    plot(lon, lat); 
end
