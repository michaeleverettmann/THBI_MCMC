%% 3d topography map
addpath('/Users/igamadan/Downloads/')
[x, y, z] = grdread2('/Users/igamadan/Dropbox (Brown)/RFcode/Data/CCP/alaska.grd');

figure;surf(x,y,z); shading interp
 load coastlines
clat=coastlat;
clon=coastlon;
hold on;
plot(coastlon,coastlat,'color',[0.3 0.3 0.3])
%axis([-144.5 -144 59.2 69.2]);
axis([-170 -140 60 70])
cmap = getPyPlot_cMap('gist_earth_r', 100);%'RdYlBu','Spectral'coolwarm_r'
colormap(cmap); c=colorbar('location', 'eastoutside');
set(gca,'CLim',[-500 2000])

view([85.9143069069985 13.8093789602347]);
%% mesh version
figure;mesh(x,y,z,'FaceColor','flat')

hold on;
plot(coastlon,coastlat,'color',[0.3 0.3 0.3])
%axis([-165 -141 65.2 65.7]); 
%axis([-162 -140 60.2 60.7]);
%axis([-144.5 -144 60.2 69.7]);
%axis([-152.5 -152 59.7 69.7])
%axis([-154.5 -154 59.2 69.7])
%axis([-170 -140 60 70]); 
% 
% axis([-154.5 -154 59.2 69.7]) % 154
% axis([-144.5 -144 60.2 69.7]) %144
% axis([-165 -141 60.2 60.5]) % 60.2
% axis([-165 -141 65.5 65.8]) % 65.2
 axis([-164 -141 68.2 68.5]) % 68.2

cmap = getPyPlot_cMap('Greys_r', 5000);%'RdYlBu','Spectral'coolwarm_r'
colormap(cmap); c=colorbar('location', 'eastoutside');
%set(gca,'CLim',[0 2900])
set(gca,'CLim',[-500 2900])
%set(gca,'ylim',[0 200])
%not this view([85.9143069069985 13.8093789602347]);%for N-S profile
%view([85.4138582983272 9.62645210883404]);
zz=zeros(length(lats),1);
hold on;
%scatter3(lons,lats,zz,50,zz,'ko','markerfacecolor','w')
%view([89.9392758239953 3.89272076184307]);% front view N-S
%view([5.24421901446621 -10.415709969789]);%for W-E profile
%view([-7.73847439527591 11.7618242767864]);

view([-6.95743060403416 3.00437693388134]); %60.2, 65.2 68.7
%view([85.6528545152766 2.51613213304686]); % -154 and -144 % May 23
set(gca,'zlim',[-8000 8000])
zticks([-8000 -6000 -4000 -2000 0])
zticklabels({'200','150','100','50','0'})
