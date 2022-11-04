%% fig setup
%if length(par.inv.datatypes)>1
%load('uni_lon.mat') % change line 112 too
%for uu=23:23%length(uni)
%load('uni_lat.mat')
close all
clear

ifgradient=0;
addpath('/Users/igamadan/Documents/Doc/functions/')
addpath('/Users/igamadan/Dropbox (Brown)/BJIcode/BayesianJointInv/ALASKA')
addpath('/Users/igamadan/Downloads/github_repo/')
load_faults
% cd '/Users/igamadan/Dropbox (Brown)/inversionResults_FEb/NewPhaseVelocities/CCP_62.20_-145.00/'
% res=('/Users/igamadan/Dropbox (Brown)/inversionResults_FEb/NewPhaseVelocities/CCP_62.20_-145.00/');
% mat = dir('*.mat'); 
% for q = 1:length(mat) 
%    load([res,mat(q).name]); 
% end
   

opts = delimitedTextImportOptions("NumVariables", 2, "Encoding", "UTF16-LE");
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["pts_lat", "pts_lon"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
tbl = readtable("/Users/igamadan/Downloads/Thesis_files/Chapter4_data/profile_pts_newphv", opts);

pts_lat = tbl.pts_lat;
pts_lon = tbl.pts_lon;
len=length(pts_lat);

clear opts tbl


load('/Users/igamadan/Dropbox (Brown)/final_model_z.mat')
%%% interpolate velocity models 
importVolcanoes
import_earthquakes
importWrangellSlab
kj=1;
%%
%init
%fig1 = figure(11); clf; set(gcf,'pos',[20 700 1400 900])


Z_toUse= final_model_z;%(0:0.1:300)';
fmLen=length(Z_toUse);
fm=zeros(fmLen,len);
Z=zeros(fmLen,len);
vslims = [3.2 5];

for ii=1:length(pts_lat)
    latP=pts_lat(ii);
    lonP=pts_lon(ii);
    formatspec='/Users/igamadan/Dropbox (Brown)/inversionResults_FEb/NewPhaseVelocities/CCP_%.2f_%i.00/final_model.mat';
    load(sprintf(formatspec,latP,lonP));  
    
%     if length(final_model.VSav)<fmLen
%         dif=abs(length(final_model.VSav)-fmLen);
%         final_model_VSav= [final_model.VSav;nan(1,dif)'];
%         fm(:,ii) = final_model_VSav;
%         
%         Z_VSav= [final_model.Z;nan(1,dif)'];
%         Z(:,ii) = Z_VSav;
%         
%     else
%         fm(:,ii) = final_model.VSav(1:fmLen);
%         Z(:,ii) = final_model.Z(1:fmLen);
%     end
%     
    
   
    
    fm(:,ii)=interp1(final_model.Z,final_model.VSav,Z_toUse,'linear');
    
    
    
    
   % ax = axes(fig1,'pos',[0.05 0.08 0.22 0.85]); hold on, 
    
    
    
%     idZ=find(final_model.Z==200);
%     Z=final_model.Z(1:idZ);
%     plot(ax,final_model.VSav(1:idZ),Z,'k','Linewidth',2)
      
%    set(ax,'ydir','reverse','fontsize',24,'ytick',[0:25:max(Z)],'ylim',[0 max(Z)],...
%         'color','none','xlim',vslims,'box','on','linewidth',2,'layer','top');
%    % if ii>1, set(ax(ii),'yticklabel',[]); end
%     %title(ax,sprintf('\\textbf{%s %s}',latP,lonP),'fontsize',25,'fontweight','bold','interpreter','latex')
%     xlabel(ax,'$\mathbf{V_S}$ \textbf{(km/s)}','Fontsize',28,'interpreter','latex')
%     
%     ylabel(ax,'\textbf{Depth (km)}','Fontsize',27,'interpreter','latex')
% 
%     lowerB =130; upperB=10;
%     lowerPVG =225; upperPVG=30;
%     %
 %   % what type of dv/dz? 1 = standard, 0 = between min and max
%    option=1; %find dv/dz
%    [Z_PVGmax,Z_PVGmin,dv_dz,pvg] = PVG_calc(final_model,lowerPVG,upperPVG,ax,1,option);
%    idx_PVGmin = find(Z==Z_PVGmin);
%    idx_PVGmax = find(Z==Z_PVGmax);
%    DvDz(ii)=dv_dz;
%    PVG_depth(ii)=pvg;
      
end

%uni=unique(pts_lon);
%%


id_toPLot=find(pts_lon==-144); 
 

%%%%%%%
lon_p=linspace(min(pts_lon),max(pts_lon),len);
lat_p=linspace(min(pts_lat),max(pts_lat),len);

[Xg, Yg] = meshgrid(lon_p, lat_p);



%%% find the profile
%  lons=Xg(id_toPLot);
%  lats=Yg(id_toPLot);
 lons=pts_lon(id_toPLot);
 lats=pts_lat(id_toPLot);


for i=1:length(event_lon)
    if event_depth(i)>50
    event_Lon(i) = event_lon(i);
    event_Lat(i) = event_lat(i);
    event_Depth(i) =event_depth(i);
    end
end

%%
figure(1)

lon=[-170 -130]; lat=[58 72];
%lat=[min(lats)-2 max(lats)+2]; lon=[min(lons)-2 max(lons)+2];%[26 40];

latlons=[min(lats) min(lons) max(lats) max(lons)];


for i=1:length(lons)
    azdiff = azimuth([latlons(kj,1),latlons(kj,2)],[latlons(kj,3),latlons(kj,4)])-...
    azimuth([latlons(kj,1),latlons(kj,2)],[lats(i),lons(i)]);
    disttmp = distance([latlons(kj,1),latlons(kj,2)],[lats(i),lons(i)]);
    dalj(i) = 111.16*atand(cosd(azdiff)*tand(disttmp));
end



daspect([111.16, 111.16*distance(mean(lat),0,mean(lat),1), 1]);
box on;  axis([lon, lat]);%[-90 -77, lat]);

%   contourf(tmplon(:,1),tmplat(1,:),log(rweighted_n_events(:,:,200))','linestyle','none');

scatter3(event_Lon, event_Lat,event_Depth,10,'k.');hold on;

set(gca,'xgrid','on','ygrid','on','layer','top')


plot(long_volcano,lat_volcano,'^','MarkerSize',10,...
'MarkerEdgeColor','red',...
'MarkerFaceColor',[1 .6 .6])
plot(lons,lats,'k-','linewidth',2); hold on;
plot(lons(1:1:end),lats(1:1:end),'ko','markerfacecolor','w','markersize',8);

load coastlines
clat=coastlat;
clon=coastlon;
view([-0.899999999999999 90])

axis([-170 -138 56 72]);
hold on
plot(coastlon,coastlat,'color',[0.3 0.3 0.3])
%colorbar('location', 'southoutside');
% text(0.85, -10.5,['\fontsize{14} \bf ' num2str(latlons(kj,1)) 'N, ' num2str(latlons(kj,2)) 'E   to   ' ...
%     num2str(latlons(kj,3)) 'N, ' num2str(latlons(kj,4)) 'W'])

 ofile=sprintf([num2str(latlons(kj,1)),'N',num2str(latlons(kj,2)),'E', ...
    num2str(latlons(kj,3)),'N',num2str(latlons(kj,4)),'W_mapLocation.jpg']);
save2jpg(1,ofile,'./');

%%
if ifgradient

figure(2)
%surf(fm(:,id_toPLot)
%surf(Xg,Yg,Zg);shading interp 
grad=gradient(fm(:,id_toPLot));
surf(dalj,Z_toUse,grad);%fm(:,id_toPLot));
shading interp %fm_line=[fm(:,2),fm(:,21),fm(:,40),fm(:,56),fm(:,57)];
cmap = getPyPlot_cMap('Spectral', 1000);%'RdYlBu','Spectral'coolwarm_r'
colormap(cmap); c=colorbar('location', 'westoutside');
xlabel(c,'Gradient','fontsize',20);
set(gca,'ylim',[0 200])
view(gca,[-180.180314158076 -90]);
set(gca,'ydir','reverse','xdir','reverse','fontsize',20,'ytick',[0:25:length(Z_toUse)],...%'ylim',[0 max(Z_toUse)],...
    'color','none','box','on','linewidth',2,'layer','top');
ylabel('Depth (km)','fontsize',20);
yticks([0 50 100 150 200 250])
yticklabels({'0','50','100','150','200','250'});hold on;
plot(dalj(1:1:end),0.*dalj(1:1:end),'ko','markersize',8,'markerfacecolor','w');

%xticks(1:length(pts_lon(id_toPLot))) %%% change here for different orientation
%xticklabels(pts_lon(id_toPLot)) %%% change here for different orientation
%view(gca,[-0.0396578538102119 90]);
%xlabel('Profile Distance (km)','fontsize',24);
 set(gca,'xlim',[0 max(dalj)]);hold on;
 
text(0.85, -10.5,['\fontsize{14} \bf ' num2str(latlons(kj,1)) 'N, ' num2str(latlons(kj,2)) 'E   to   ' ...
    num2str(latlons(kj,3)) 'N, ' num2str(latlons(kj,4)) 'W'])

ofile2=sprintf([num2str(latlons(kj,1)),'N',num2str(latlons(kj,2)),'E', ...
    num2str(latlons(kj,3)),'N',num2str(latlons(kj,4)),'gradient.jpg']);

for i=1:length(lon_f) 
    azdiff = azimuth([latlons(kj,1),latlons(kj,2)],[latlons(kj,3),latlons(kj,4)])-...
        azimuth([latlons(kj,1),latlons(kj,2)],[lat_f(i),lon_f(i)]);
    disttmp = distance([latlons(kj,1),latlons(kj,2)],[lat_f(i),lon_f(i)]);
    disteq(i) = 111.16*atand(cosd(azdiff)*tand(disttmp));
    disteq(i) = disteq(i);
    distnorm(i) = 111.16*atand(sind(azdiff)*tand(disttmp));       
        if abs(distnorm(i)) < 2  %distance from profiles to EQ in km
            plot (disteq(i),0,'k|','LineWidth',3,'markersize',20); hold on;
        end          
end  
save2jpg(2,ofile2,'./');
end
 %%
%% SURF PLOT

Zg=imgaussfilt(fm(:,id_toPLot),0.7);

figure
%surf(fm(:,id_toPLot)
%surf(Xg,Yg,Zg);shading interp 
surf(dalj,Z_toUse,Zg);
%surf(dalj,Z_toUse,Zg);
shading interp %fm_line=[fm(:,2),fm(:,21),fm(:,40),fm(:,56),fm(:,57)];
cmap = getPyPlot_cMap('Spectral', 1000);%'RdYlBu','Spectral'coolwarm_r'
colormap(cmap); c=colorbar('location', 'westoutside');
xlabel(c,'Average Vs Velocity (km/s)','fontsize',20);
%set(gca,'ylim',[0 250],'CLim',[3.2 4.7])
set(gca,'ylim',[0 200],'CLim',[4.1 4.66])
view(gca,[-180.180314158076 -90]);
set(gca,'ydir','reverse','xdir','reverse','fontsize',20,'ytick',[0:25:length(Z_toUse)],...%'ylim',[0 max(Z_toUse)],...
    'color','none','box','on','linewidth',2,'layer','top');
ylabel('Depth (km)','fontsize',20);
yticks([0 50 100 150 200 250])
yticklabels({'0','50','100','150','200','250'});hold on;
plot(dalj(1:1:end),0.*dalj(1:1:end),'ko','markersize',8,'markerfacecolor','w');

%xticks(1:length(pts_lon(id_toPLot))) %%% change here for different orientation
%xticklabels(pts_lon(id_toPLot)) %%% change here for different orientation
%view(gca,[-0.0396578538102119 90]);
%xlabel('Profile Distance (km)','fontsize',24);
 set(gca,'xlim',[0 max(dalj)]);hold on;
 
text(0.85, -10.5,['\fontsize{14} \bf ' num2str(latlons(kj,1)) 'N, ' num2str(latlons(kj,2)) 'E   to   ' ...
    num2str(latlons(kj,3)) 'N, ' num2str(latlons(kj,4)) 'W'])

 ofile3=sprintf([num2str(latlons(kj,1)),'N',num2str(latlons(kj,2)),'E', ...
    num2str(latlons(kj,3)),'N',num2str(latlons(kj,4)),'W_zoomIn.jpg']);
%save2jpg(3,ofile,'./');
%
 

 
  for i=1:length(lat_volcano)
      azdiff_vol = azimuth([latlons(kj,1),latlons(kj,2)],[latlons(kj,3),latlons(kj,4)])-...
            azimuth([latlons(kj,1),latlons(kj,2)],[lat_volcano(i),long_volcano(i)]);
        disttmp_vol = distance([latlons(kj,1),latlons(kj,2)],[lat_volcano(i),long_volcano(i)]);
        distvol(i) = 111.16*atand(cosd(azdiff_vol)*tand(disttmp_vol));
        distvol(i) = distvol(i);
        distnorm_vol(i) = 111.16*atand(sind(azdiff_vol)*tand(disttmp_vol));       
        if abs(distnorm_vol(i)) < 50  %distance from profiles to EQ in km
            plot (distvol(i),0,'^','MarkerSize',30,...
            'MarkerEdgeColor','red',...
            'MarkerFaceColor',[1 .6 .6]); hold on;
        end   
        
  end
   

  %

for i=1:length(event_Lon) 
    azdiff = azimuth([latlons(kj,1),latlons(kj,2)],[latlons(kj,3),latlons(kj,4)])-...
        azimuth([latlons(kj,1),latlons(kj,2)],[event_Lat(i),event_Lon(i)]);
    disttmp = distance([latlons(kj,1),latlons(kj,2)],[event_Lat(i),event_Lon(i)]);
    disteq(i) = 111.16*atand(cosd(azdiff)*tand(disttmp));
    disteq(i) = disteq(i);
    distnorm(i) = 111.16*atand(sind(azdiff)*tand(disttmp));       
        if abs(distnorm(i)) < 50  %distance from profiles to EQ in km
            plot (disteq(i),event_Depth(i),'ko','markerfacecolor','w','markersize',3); hold on;
        end          
end  

for i=1:length(lon_f) 
    azdiff = azimuth([latlons(kj,1),latlons(kj,2)],[latlons(kj,3),latlons(kj,4)])-...
        azimuth([latlons(kj,1),latlons(kj,2)],[lat_f(i),lon_f(i)]);
    disttmp = distance([latlons(kj,1),latlons(kj,2)],[lat_f(i),lon_f(i)]);
    disteq(i) = 111.16*atand(cosd(azdiff)*tand(disttmp));
    disteq(i) = disteq(i);
    distnorm(i) = 111.16*atand(sind(azdiff)*tand(disttmp));       
        if abs(distnorm(i)) < 2  %distance from profiles to EQ in km
            plot (disteq(i),0,'k|','LineWidth',3,'markersize',20); hold on;
        end          
end  
%plot(lon_f,lat_f,'k','LineWidth',2)




 for i=1:length(WraLon) 
    azdiff = azimuth([latlons(kj,1),latlons(kj,2)],[latlons(kj,3),latlons(kj,4)])-...
        azimuth([latlons(kj,1),latlons(kj,2)],[WraLat(i),WraLon(i)]);
    disttmp = distance([latlons(kj,1),latlons(kj,2)],[WraLat(i),WraLon(i)]);
    disteqWra(i) = 111.16*atand(cosd(azdiff)*tand(disttmp));
    disteqWra(i) = disteqWra(i);
    distnormWra(i) = 111.16*atand(sind(azdiff)*tand(disttmp));       
  %  plot(disteqWra,WraDepth(i),'k.','markerfacecolor','w','markersize',10); hold on;
    if abs(distnormWra(i)) < 20  %distance from profiles to EQ in km
        plot(disteqWra(i),WraDepth(i),'ko','markerfacecolor','w','markersize',4); hold on;
    end          
 end  
save2jpg(3,ofile3,'./');
%clearvars -except uu uni
%close all

% 