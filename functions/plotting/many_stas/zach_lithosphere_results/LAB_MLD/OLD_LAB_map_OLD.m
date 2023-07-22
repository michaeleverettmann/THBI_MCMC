clear all
%% load in 3D grid
load('surface_colated_b1_V7.mat');
% establish grid size etc.
[Nx,Ny] = size(xgrid);
Nz = length(depths);

smthnpts = 5;
smthdist = 50;
ifsave = false;

mldmax = 120;
minmlddv = 0.004;

% %% test by plotting moho
figure(901); clf
surface(xgrid,ygrid,zmoh_surf)
colorbar,colormap(flipud(parula))

addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')


%% loop over all the vertical profiles and find lab/mld
zlab_surf_stat = nan(size(zmoh_surf));
zlab_surf_vmi = nan(size(zmoh_surf));
zlab_surf_v50 = nan(size(zmoh_surf));
zlab_surf_v90 = nan(size(zmoh_surf));
zmld_surf_stat = nan(size(zmoh_surf));
zmld_surf_vmi = nan(size(zmoh_surf));

%     metha = ["stat","vmi",]

% zlab_surf_all

for ix = 1:Nx
for iy = 1:Ny
    Vs = squeeze(mgrid3d(ix,iy,:));
    try

        [lab_provisional,mld_provisional] = ...
            LAB_finder(Vs,depths,zmoh_surf(ix,iy),smthnpts,0,minmlddv);
   catch
        try
            [lab_provisional,mld_provisional] = ...
            LAB_finder(Vs,depths,zmoh_surf(ix,iy),smthnpts,0,280);
        catch
            fprintf('Some failure for [%.2f N, %.2f E]\n',latgrid(ix,iy),longrid(ix,iy))
            lab_provisional = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
            mld_provisional = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
        end
    end
%     continue
    % issue with mld mistaken for lab
    if mld_provisional.z_stat > mldmax
        lab_provisional = mld_provisional; 
        mld_provisional = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
    end
    metha = ["stat","vmi",]

    % now insert values if not empty/nan
    if ~isempty(lab_provisional.z_stat) && ~isnan(lab_provisional.z_stat)
        zlab_surf_stat(ix,iy) = lab_provisional.z_stat;
    end
    if ~isempty(lab_provisional.z_vmi) && ~isnan(lab_provisional.z_vmi)
        zlab_surf_vmi(ix,iy) = lab_provisional.z_vmi;
    end
    if ~isempty(lab_provisional.z_vmi50ma) && ~isnan(lab_provisional.z_vmi50ma)
        zlab_surf_v50(ix,iy) = lab_provisional.z_vmi50ma;
    end
    if ~isempty(lab_provisional.z_vmi90ma) && ~isnan(lab_provisional.z_vmi90ma)
        zlab_surf_v90(ix,iy) = lab_provisional.z_vmi90ma;
    end
%     if mld_provisional.z_stat > 100
%         keyboard, 
%     end
    if ~isempty(mld_provisional.z_stat) && ~isnan(mld_provisional.z_stat)
        zmld_surf_stat(ix,iy) = mld_provisional.z_stat;  
    end
    if ~isempty(mld_provisional.z_vmi) && ~isnan(mld_provisional.z_vmi)
        zmld_surf_vmi(ix,iy) = mld_provisional.z_vmi;  
    end

end
end

% kill points too far from real stations - not robust
load("_distance_pt_to_sta.mat");
zlab_surf_stat(pt_dist_km>70) = nan;
zlab_surf_vmi(pt_dist_km>70) = nan;
zlab_surf_v50(pt_dist_km>70) = nan;
zlab_surf_v90(pt_dist_km>70) = nan;
zmld_surf_stat(pt_dist_km>70) = nan;
zmld_surf_vmi(pt_dist_km>70) = nan;
zmoh_surf(pt_dist_km>70) = nan;

%% LAB establish mean depths, and some simple stats
fprintf('\nLAB info\n')
meth = ["stat","v50","v90","vmi"];
zlab_mean = zeros(4,1);
for im = 1:length(meth)
    fprintf('%s method:\n',meth(im))
    eval(['zlab_mean(im) = nanmean(zlab_surf_',char(meth(im)),',''all'');'])
    fprintf('  Mean depth = %.1f km\n',zlab_mean(im))
    fprintf('  Npoints = %.0f\n',eval(['sum(~isnan(zlab_surf_',char(meth(im)),'),''all'')']))
end


%concatenate LABs
zlab_surf_all = cat(3,zlab_surf_stat,zlab_surf_v50,zlab_surf_v90,zlab_surf_vmi);
% get overall mean LAB

% get differential LAB depths
zlab_diff_all = cat(3,zlab_surf_stat - nanmean(zlab_surf_stat,'all'),...
                      zlab_surf_v50 - nanmean(zlab_surf_v50,'all'),...
                      zlab_surf_v90 - nanmean(zlab_surf_v90,'all'),...
                      zlab_surf_vmi - nanmean(zlab_surf_vmi,'all'));

% make map of lab uncertainty = proportional to how much zlab estimates
% diverge
zlab_std_all = std(zlab_surf_all,1,3);

% decide on mask - nan out if more than one method gives nan lab
zlab_nan = sum(isnan(zlab_surf_all),3)>1;

% plot uncertainty
figure(76);clf
pcolor(longrid,latgrid,zlab_std_all)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45])
colorbar,colormap(flipud(parula(15))),caxis([30 50])

% get smoothed absolute LAB depths
zlab_surfsmooth_all = nan(size(zlab_surf_all));
for im = 1:4
    for ix = 1:Nx
    for iy = 1:Ny
        if pt_dist_km(ix,iy)>70, continue; end
        zlab_surfsmooth_all(ix,iy,im) = ...
            gaussian_smoothing_lab(ix,iy,longrid,latgrid,zlab_surf_all(:,:,im),1./zlab_std_all,smthdist);
    end
    end
end

% get smoothed differential LAB depths
zlab_diffsmooth_all = nan(size(zlab_diff_all));
for im = 1:4
    for ix = 1:Nx
    for iy = 1:Ny
        if pt_dist_km(ix,iy)>70, continue; end
        zlab_diffsmooth_all(ix,iy,im) = ...
            gaussian_smoothing_lab(ix,iy,longrid,latgrid,zlab_diff_all(:,:,im),1./zlab_std_all,smthdist);
    end
    end
end

% Make preferred map - use average differential value from all the methods,
% pin to mean of v50, as a halfway house

zlab_pref = zlab_mean(2) + mean(zlab_diffsmooth_all,3);

%% MLD establish mean depths, and some simple stats
fprintf('\nMLD info\n')
methm = ["stat","v50","v90","vmi"];
zmld_mean = zeros(length(methm),1);
for im = 1:length(methm)
    fprintf('%s method:\n',methm(im))
    eval(['zmld_mean(im) = nanmean(zmld_surf_',char(methm(im)),',''all'');'])
    fprintf('  Mean depth = %.1f km\n',zmld_mean(im))
    fprintf('  Npoints = %.0f\n',eval(['sum(~isnan(zmld_surf_',char(methm(im)),'),''all'')']))
end


%concatenate mlds
zmld_surf_all = cat(3,zmld_surf_stat,zmld_surf_vmi);
% get overall mean mld

% get differential mld depths
zmld_diff_all = cat(3,zmld_surf_stat - nanmean(zmld_surf_stat,'all'),...
                      zmld_surf_vmi - nanmean(zmld_surf_vmi,'all'));

% make map of mld uncertainty = proportional to how much zmld estimates
% diverge
zmld_std_all = abs(zmld_surf_stat - zmld_surf_vmi);

% decide on mask - nan out if even one method gives nan mld
zmld_nan = sum(isnan(zmld_surf_all),3)>0;

% plot uncertainty
figure(76);clf
pcolor(longrid,latgrid,zmld_std_all)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45])
colorbar,colormap(flipud(parula(15))),caxis([30 50])

% get smoothed absolute mld depths
zmld_surfsmooth_all = nan(size(zmld_surf_all));
for im = 1:2
    for ix = 1:Nx
    for iy = 1:Ny
        if pt_dist_km(ix,iy)>70, continue; end
        zmld_surfsmooth_all(ix,iy,im) = ...
            gaussian_smoothing_lab(ix,iy,longrid,latgrid,zmld_surf_all(:,:,im),1./zmld_std_all,smthdist);
    end
    end
end

% get smoothed differential mld depths
zmld_diffsmooth_all = nan(size(zmld_diff_all));
for im = 1:2
    for ix = 1:Nx
    for iy = 1:Ny
        if pt_dist_km(ix,iy)>70, continue; end
        zmld_diffsmooth_all(ix,iy,im) = ...
            gaussian_smoothing_lab(ix,iy,longrid,latgrid,zmld_diff_all(:,:,im),1./zmld_std_all,smthdist);
    end
    end
end

% Make preferred map - use average differential value from all the methods,
% pin to mean of v50, as a halfway house

zmld_pref = mean(zmld_surfsmooth_all,3);


%% PLOTTING
addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')



%moho
figure(901); clf, hold on
surf(longrid,latgrid,zmoh_surf)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45])
colorbar,colormap(flipud(parula(15))),caxis([30 50])

%lab
figure(902); clf, hold on
surface(longrid,latgrid,zlab_surf_stat)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45]), title('LAB - zstat')
colorbar,colormap(flipud(turbo(15))),caxis([60 280])

%mld
figure(903); clf, hold on
surface(longrid,latgrid,zmld_surf_stat)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45]),title('MLD - zstat')
colorbar,colormap(flipud(turbo(15))),caxis([50 120])

%lab2
figure(912); clf, hold on
surface(longrid,latgrid,zlab_surf_vmi)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45]), title('LAB - vmi')
colorbar,colormap(flipud(turbo(15))),caxis([70 280])

%mld2
figure(913); clf, hold on
surface(longrid,latgrid,zmld_surf_vmi)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45]), title('MLD - vmi')
colorbar,colormap(flipud(turbo(15))),caxis([50 120])



% un-smoothed diff-LAB maps, 4 methods
figure(444),clf
for im = 1:4
    subplot(2,2,im);
    surf(longrid,latgrid,zeros(size(longrid)),zlab_diff_all(:,:,im),'linestyle','none')
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord)
        line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
    end
    xlim([-90 -70]);ylim([27,45]), view(0,90)
    colorbar,colormap(flipud(parula(15))),caxis([-30 30])
    title(meth(im),'fontsize',30);
end

% SMOOTHed diff-LAB maps, 4 methods
figure(445),clf
for im = 1:4
    subplot(2,2,im); hold on
    ZZ = zlab_diffsmooth_all(:,:,im);
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.2)
    ZZ(isnan(zlab_surf_all(:,:,im))) = nan;
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord)
        line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
    end
    xlim([-90 -70]);ylim([27,45]), view(0,90)
    colorbar,colormap(flipud(parula(15))),caxis([-40 40])
    title(meth(im),'fontsize',30);
end

% look at histograms for range of differential values
figure(446),clf
for im = 1:4
    subplot(2,2,im);
    aa = zlab_diffsmooth_all(:,:,im); aa = aa(:);
    hist(aa)
  
    title(meth(im),'fontsize',30);
end

%% LAB nice figure
% three columns, two rows
%  - left shows min LAB (stat) smoothed and differential
%  - right shows max LAB (vmi) smoothed and differential
%  - middle shows pref LAB (mean-v50 + mean(diff_lab)) smoothed and differential
% axpref(row,col)

figure(489),set(gcf,'position',[512 465 1050 872]),clf
axw = 0.265; axh = 0.4;
axpref(1,1) = axes('position',[0.04 0.52 axw axh]); hold on
axpref(1,2) = axes('position',[0.34 0.52 axw axh]); hold on
axpref(1,3) = axes('position',[0.64 0.52 axw axh]); hold on
axpref(2,1) = axes('position',[0.04 0.07 axw axh]); hold on
axpref(2,2) = axes('position',[0.34 0.07 axw axh]); hold on
axpref(2,3) = axes('position',[0.64 0.07 axw axh]); hold on

% STAT-ABS
ZZ = zlab_surfsmooth_all(:,:,1);
surf(axpref(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zlab_surf_all(:,:,1))) = nan;
surf(axpref(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% PREF-ABS
ZZ = zlab_pref;
surf(axpref(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(zlab_nan) = nan;
surf(axpref(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% VMI-ABS
ZZ = zlab_surfsmooth_all(:,:,4);
surf(axpref(1,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zlab_surf_all(:,:,4))) = nan;
surf(axpref(1,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% STAT-DIFF
ZZ = zlab_diffsmooth_all(:,:,1);
surf(axpref(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zlab_surf_all(:,:,1))) = nan;
surf(axpref(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% PREF-DIFF
ZZ = zlab_pref-zlab_mean(2);
surf(axpref(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(zlab_nan) = nan;
surf(axpref(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% VMI-DIFF
ZZ = zlab_diffsmooth_all(:,:,4);
surf(axpref(2,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zlab_surf_all(:,:,4))) = nan;
surf(axpref(2,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% adjust colourmaps
for ia = 1:3
    colormap(axpref(1,ia),flipud(parula(20))),caxis(axpref(1,ia),[90 250])
    colormap(axpref(2,ia),flipud(turbo(20))),caxis(axpref(2,ia),[-60 60])
end
% put on geography
for imp = 1:numel(axpref)    
    geogax(axpref(imp));
end
% put on CAMP
a = load('../CAMP_digitize/CAMP_fromGao2000f1_nan.mat');
CAMP = a.a;
for imp = 1:numel(axpref)    
    plot(axpref(imp),CAMP(:,2),CAMP(:,3),'w','Linewidth',1.5);
end

%ticks
set(axpref(:,[2,3]),'yticklabel','')
% titles
titles = {'LAB depth (min)','LAB depth (pref)','LAB depth (max)'};
for ia = 1:3
    title(axpref(1,ia),titles{ia},'fontsize',25,'fontweight','bold')
end
% colorbars
cbtitles = {'Absolute depth (km)','Relative depth (km)'};
for ii = 1:2
    hcb(ii) = colorbar(axpref(ii,3),'eastoutside','fontsize',15); drawnow
    ap = axpos(axpref(ii,3)); ap(3) = axpos(axpref(ii,2),3);
    set(axpref(ii,3),'position',ap)
    ylabel(hcb(ii),cbtitles{ii},'fontsize',18,'fontweight','bold')
end
set(gcf,'color','w')
return
if ifsave
    save2pdf(489,'LAB_maps');
end


%% MLD nice figure
% two columns, two rows
%  - left shows stat MLD smoothed and differential
%  - right shows vmi MLD smoothed and differential
% axpref(row,col)

figure(490),set(gcf,'position',[1355 465 740 872]),clf
axw = 0.39; axh = 0.4;
axprefm(1,1) = axes('position',[0.04 0.52 axw axh]); hold on
axprefm(1,2) = axes('position',[0.48 0.52 axw axh]); hold on
axprefm(2,1) = axes('position',[0.04 0.07 axw axh]); hold on
axprefm(2,2) = axes('position',[0.48 0.07 axw axh]); hold on

% STAT-ABS
ZZ = zmld_surfsmooth_all(:,:,1);
% surf(axprefm(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zmld_surf_all(:,:,1))) = nan;
surf(axprefm(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% VMI-ABS
ZZ = zmld_surfsmooth_all(:,:,2);
% surf(axprefm(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zmld_surf_all(:,:,2))) = nan;
surf(axprefm(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% STAT-DIFF
ZZ = zmld_diffsmooth_all(:,:,1);
% surf(axprefm(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zmld_surf_all(:,:,1))) = nan;
surf(axprefm(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% VMI-DIFF
ZZ = zmld_diffsmooth_all(:,:,2);
% surf(axprefm(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(zmld_surf_all(:,:,2))) = nan;
surf(axprefm(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% adjust colourmaps
for ia = 1:2
    colormap(axprefm(1,ia),flipud(parula(20))),caxis(axprefm(1,ia),[50 120])
    colormap(axprefm(2,ia),flipud(turbo(20))),caxis(axprefm(2,ia),[-20 20])
end
% put on geobgraphy
for imp = 1:numel(axprefm)    
    geogax(axprefm(imp));
end
%ticks
set(axprefm(:,2),'yticklabel','')
% titles
titles = {'MLD depth (stat)','MLD depth (vmi)'};
for ia = 1:2
    title(axprefm(1,ia),titles{ia},'fontsize',25,'fontweight','bold')
end
% colorbars
cbtitles = {'Absolute depth (km)','Relative depth (km)'};
for ii = 1:2
    hcb(ii) = colorbar(axprefm(ii,2),'eastoutside','fontsize',15); drawnow
    ap = axpos(axprefm(ii,2)); ap(3) = axpos(axprefm(ii,1),3);
    set(axprefm(ii,2),'position',ap)
    ylabel(hcb(ii),cbtitles{ii},'fontsize',18,'fontweight','bold')
end
set(gcf,'color','w')

if ifsave
    save2pdf(490,'MLD_maps');
end


%%
for im = 1:4
    subplot(2,2,im); hold on
    ZZ = zlab_diffsmooth_all(:,:,im);
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.2)
    ZZ(isnan(zlab_surf_all(:,:,im))) = nan;
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord)
        line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
    end
    xlim([-90 -70]);ylim([27,45]), view(0,90)
    colorbar,colormap(flipud(parula(15))),caxis([-40 40])
    title(meth(im),'fontsize',30);
end

    

%% saving
LAB_info = struct('longrid',longrid,'latgrid',latgrid,'gauss_smth_dist_km',smthdist,'Vs_smthpts',smthnpts,...
                  'zlab_pref',zlab_pref,'zdlab_pref',mean(zlab_diffsmooth_all,3),...
                  'zlab_isnan',zlab_nan,'zlab_wt',1./zlab_std_all,...
                  'zlab_surfsmooth_all',zlab_surfsmooth_all,'zlab_methods',meth);


MLD_info = struct('longrid',longrid,'latgrid',latgrid,'gauss_smth_dist_km',smthdist,'Vs_smthpts',smthnpts,...
                  'zmld_surf_all',zmld_surf_all,'zlab_methods',methm);

save('LAB_result','LAB_info','-v7.3')
save('MLD_result','MLD_info','-v7.3')

%% do a little post-processing
smthnpts = 20;
for ix = 2:Nx
for iy = 1:Ny
    if isnan(zlab_surf_stat(ix,iy)), continue; end
    if abs(zlab_surf_stat(ix,iy)-zlab_surf_stat(ix-1,iy)) > 20
        
        Vs1 = squeeze(mgrid3d(ix,iy,:));
        Vs0 = squeeze(mgrid3d(ix-1,iy,:));

        figure(902)
        hpt = plot(longrid(ix,iy),latgrid(ix,iy),'p');
        
        figure(434),clf, hold on
        plot(smooth(depths,Vs1,smthnpts),depths(:),'r','linewidth',2);
        plot(smooth(depths,Vs0,smthnpts),depths(:),'b','linewidth',2);
        legend('Vs1','Vs0')
        set(gca,'ydir','reverse'); hold on
        yline(zlab_surf_stat(ix,iy),'--r','linewidth',2)
        yline(zlab_surf_stat(ix-1,iy),'--b','linewidth',2)


        [lab_provisional,mld_provisional] = ...
            LAB_finder(Vs1,depths,zmoh_surf(ix,iy),smthnpts,1);
        clone_figure(99,100)
        [lab_provisional,mld_provisional] = ...
            LAB_finder(Vs0,depths,zmoh_surf(ix,iy),smthnpts,1);
    end
end
end


%% make some smoothed maps
zlab_surfsmth_stat = nan(size(zlab_surf_stat));
for ix = 1:Nx
for iy = 1:Ny
    zlab_surfsmth_stat(ix,iy) = gaussian_smoothing_lab(ix,iy,longrid,latgrid,zlab_surf_stat,1./zlab_std_all,50);
end
end
zlab_surfsmth_stat(pt_dist_km>70) = nan;


figure(454);clf

subplot(1,2,1);
surf(longrid,latgrid,zlab_surf_stat)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45]), view(0,90)
colorbar,colormap(flipud(parula(15)))
title('RAW');

subplot(1,2,2);
surf(longrid,latgrid,zlab_surfsmth_stat)
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim([-90 -70]);ylim([27,45]), view(0,90)
colorbar,colormap(flipud(parula(15)))
title('SMOOTHED');




%% subfunctions
function z_smth = gaussian_smoothing_lab(ix,iy,longrid,latgrid,zgrid,wtgrid,smthdist)
% function to apply a gaussian smoothing function across a surface,
% smoothing point ix, iy (i.e. grid(ix,iy)) using the points around it.
% smthdist is the sigma for the gaussian (i.e. 1/4 of the 95% interval

d2k = 111.1949;
% find distance IN KM to all other points and use to weight
dX = abs(longrid - longrid(ix,iy))*d2k*cosd(latgrid(ix,iy)); % account for sphericity in deg2km
dY = abs(latgrid - latgrid(ix,iy))*d2k;
dR2 = dX.^2 + dY.^2; % square distance, in km^2
dRwt = exp(-dR2./(2*smthdist.^2));

nnan = ~isnan(zgrid) & ~isnan(wtgrid);

% calculate the smoothed value
z_smth = sum(zgrid(nnan).*dRwt(nnan).*wtgrid(nnan))./sum(dRwt(nnan).*wtgrid(nnan));
end

function ax = geogax(ax)
addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord)
    line(ax,lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
end
xlim(ax,[-90 -70]);
ylim(ax,[27,45]);
view(ax,0,90)
set(ax,'fontsize',15,'linewidth',2,'box','on','layer','top','Xgrid','off','Ygrid','off')
end


