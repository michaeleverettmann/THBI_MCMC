clear 
%% load all important results
% fast lith
load('../LAB_MLD/fastlith_map.mat')
% LAB from temperature analysis
T_lab = 1150;
load(['../LAB_MLD/LAB_T',num2str(T_lab),'.mat'])
% MLD from gradient analysis
load('../LAB_MLD/MLD_vgrads.mat')
% geologic fronts
appf = appalachian_front;
gref = grenville_front;

ifsave = true;

%% necessary paths
addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
% addpath('/Users/zeilon/Dropbox/MATLAB/lib/m_map/')
addpath('~/MATLAB/m_map');
addpath('/Users/brennanbrunsvik/Documents/repositories/Base_code/colormaps/colorbrewer')

%% Make figure with main lithospheric speed/thickness results
figure(1),clf,set(gcf,'pos',[40 250 1168 1081],'color','w')
ax11 = axes('pos',[0.0558 0.5109 0.4 0.44]);hold on;m_geogax(ax11);
ax12 = axes('pos',[0.54 0.5109 0.4 0.44]);hold on;m_geogax(ax12);
ax21 = axes('pos',[0.0558 0.0354 0.4 0.44]);hold on;m_geogax(ax21);
ax22 = axes('pos',[0.54 0.0354 0.4 0.44]);hold on;m_geogax(ax22);

%% plot T-based LAB
axes(ax11)
m_contourf(longrid,latgrid,z_lab_Tiso_smth,30,'linestyle','none')
m_geogax(ax11);

title(['Z_{LAB} (km) for T = ',num2str(T_lab),' C'])
colormap(ax11,brewermap(15,'PuOr'))
caxis(ax11,[80 260])


% % % % put on CAMP
% % % a = load('../CAMP_digitize/CAMP_fromGao2000f1_nan.mat');
% % % CAMP = a.a;
% % % m_plot(CAMP(:,2),CAMP(:,3),'k','Linewidth',2);
% % % m_plot(CAMP(:,2),CAMP(:,3),'g','Linewidth',1);

% fronts
m_plot(appf(:,1),appf(:,2),'--','linewidth',2,'color',[44,162,95]/255)
m_plot(gref(:,1),gref(:,2),'--','linewidth',2,'color',[153,206,201]/255)

% add colourbar
ax11_ = axes('Position',axpos(ax11),'visible','off','xlim',[0 1],'ylim',[0 1]);
cbar_custom(ax11_,'location',[0.8 0.85 0.1 0.6],'lims',get(ax11,'CLim'),...
            'cmap',get(ax11,'Colormap'),'ncols',size(get(ax11,'Colormap'),1),...
            'tickside','right','tickvals',[100:25:300],'fontsize',14,...
            'title','LAB depth (km)')

%% plot gradient-based MLD
% calculate size scale on basis of depth
zz = 51:150; % depth
sz = 1000*log10(zz./50).^2; % associated size

axes(ax12);
m_scatter(longrid(:),latgrid(:),interp1(zz,sz,MLD_info.zmld_pref(:)),MLD_info.vmld_pref(:),...
    'filled','markeredgecolor',0.2*[1 1 1],'markerfacealpha',1)
title('MLD - pref')
colormap(ax12,brewermap(15,'BrBG')),caxis(ax12,[4.3 4.65]);

% fronts
m_plot(appf(:,1),appf(:,2),'--','linewidth',2,'color',[197,27,138]/255)
m_plot(gref(:,1),gref(:,2),'--','linewidth',2,'color',[221,28,119]/255)
% add scale for size
axes(ax12)
zz_scale = [60:10:110];
lon_scale = -72.3;
for iz = 1:length(zz_scale)
    m_scatter(lon_scale,39-iz,interp1(zz,sz,zz_scale(iz)),1,'markeredgecolor',0.2*[1 1 1],'markerfacecolor',0.5*[1 1 1]);
    m_text(lon_scale-0.7,39-iz,[num2str(zz_scale(iz)),' km'],'horizontalalignment','right','verticalalignment','middle','fontsize',12)
end
% add colourbar
ax12_ = axes('Position',axpos(ax12),'visible','off','xlim',[0 1],'ylim',[0 1]);
cbar_custom(ax12_,'location',[0.87 0.92 0.1 0.6],'lims',get(ax12,'CLim'),...
            'cmap',get(ax12,'Colormap'),'ncols',size(get(ax12,'Colormap'),1),...
            'tickside','right','tickvals',[4.3:0.1:4.7],'fontsize',14,...
            'title','MLD velocity (km/s)')

%% plot integral of fast V
axes(ax21)
m_contourf(longrid,latgrid,fastlithsmth_norm,30,'linestyle','none')
% put a contour on
m_contour(longrid,latgrid,fastlithsmth_norm,1*[1 1],'linewidth',1.5,'color','k')
m_contour(longrid,latgrid,fastlithsmth_norm,1.5*[1 1],'linewidth',2,'color','k')
m_contour(longrid,latgrid,fastlithsmth_norm,2*[1 1],'linewidth',2.5,'color','k')

m_geogax(ax21);
colormap(ax21,flipud(turbo(15))),
caxis(ax21,[0.1,2])
title('FastLith (relative to cont-LProt)')

% add colourbar
ax21_ = axes('Position',axpos(ax21),'visible','off','xlim',[0 1],'ylim',[0 1]);
cbar_custom(ax21_,'location',[0.8 0.85 0.1 0.6],'lims',get(ax21,'CLim'),...
            'cmap',get(ax21,'Colormap'),'ncols',size(get(ax21,'Colormap'),1),...
            'tickside','right','tickvals',[0:0.5:2],'fontsize',14,...
            'title','Relative lithosphere speed')

%% plot average UM V
axes(ax22)
m_contourf(longrid,latgrid,avvlithsmth,30,'linestyle','none')
% put a contour on
m_contour(longrid,latgrid,avvlithsmth,4.55*[1 1],'linewidth',1.5,'color','k')
m_contour(longrid,latgrid,avvlithsmth,4.6*[1 1],'linewidth',2,'color','k')
m_contour(longrid,latgrid,avvlithsmth,4.65*[1 1],'linewidth',2.5,'color','k')

m_geogax(ax22);

title('Mean UM Vs (km/s)')
colormap(ax22,flipud(parula(15))), 
caxis(ax22,[4.52 4.7])

% add colourbar
ax22_ = axes('Position',axpos(ax22),'visible','off','xlim',[0 1],'ylim',[0 1]);
cbar_custom(ax22_,'location',[0.8 0.85 0.1 0.6],'lims',get(ax22,'CLim'),...
            'cmap',get(ax22,'Colormap'),'ncols',size(get(ax22,'Colormap'),1),...
            'tickside','right','tickvals',[4.5:0.05:4.7],'fontsize',14,...
            'title','Mean UM velocity (km/s)')

%% save
if ifsave
    save2jpg(1,'Lith_analysis')
end

