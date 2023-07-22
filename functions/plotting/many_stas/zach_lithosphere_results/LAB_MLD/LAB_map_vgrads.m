clear all
%% load in 3D grid
load('surface_colated_b1_V7.mat');
% establish grid size etc.
[Nx,Ny] = size(xgrid);
Nz = length(depths);

smthnpts = 5;
smthdist = 30;
mindist2sta = 70;

ifsave = true;

mldmax = 110;
minmlddv = 0.006;
minlabdv = 0.006;

addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
addpath('~/Dropbox/MATLAB/lib/brewermap/')
addpath('../model_process_plot/')

load("_distance_pt_to_sta.mat");
zmoh_surf(pt_dist_km>mindist2sta) = nan;

%% loop over all the vertical profiles and find lab/mld
metha = ["stat","vmi50ma","vmi90ma","vmi"];
methb = ["stat","v50","v90","vmi"];
nvgs = ["lab","mld"];

% set up results structures
% depth
zlab_surf_all = nan(Nx,Ny,length(metha));
zmld_surf_all = zlab_surf_all;             
% velocity
vlab_surf_all = zlab_surf_all;
vmld_surf_all = zlab_surf_all;             
% differential velocity from max above
dvlab_surf_all = zlab_surf_all;
dvmld_surf_all = zlab_surf_all;             


for ix = 1:Nx
for iy = 1:Ny
    % don't bother if too far from a seismic station
    if pt_dist_km(ix,iy) > mindist2sta, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    % get Vs profile here.
    Vs = squeeze(mgrid3d(ix,iy,:));
    try
        [lab_prov,mld_prov] = ...
            LAB_finder(Vs,depths,zmoh_surf(ix,iy),smthnpts,0,minmlddv,minlabdv);
   catch
        try
            [lab_prov,mld_prov] = ...
            LAB_finder(Vs,depths,zmoh_surf(ix,iy),smthnpts,1,minmlddv,minlabdv);
        catch
            fprintf('Some failure for [%.2f N, %.2f E]\n',latgrid(ix,iy),longrid(ix,iy))
            lab_prov = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
            mld_prov = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
        end
    end

    % issue with mld mistaken for lab
    if mld_prov.z_stat > mldmax || mld_prov.z_vmi > mldmax+20
        lab_prov = mld_prov; 
        mld_prov = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
    end
    
    % insert into output structure
    for ityp = ["z","v","dv"]
    for invg = ["lab","mld"]
    for im = 1:length(metha)
        if ~isnan(eval(sprintf('%s_prov.z_%s',invg,metha(im))))
%         fprintf('$s%s_surf_all(ix,iy,im) = %s_prov.z_%s;\n',invg,invg,metha(im))
        eval(sprintf('%s%s_surf_all(ix,iy,im) = %s_prov.%s_%s;',ityp,invg,invg,ityp,metha(im)))
        end
    end 
    end
    end
    
%     % look at any weird ones
%     if any(isnan(zlab_surf_all(ix,iy,:)))
%         [lab_prov,mld_prov] = ...
%             LAB_finder(Vs,depths,zmoh_surf(ix,iy),smthnpts,1,minmlddv,minlabdv);
%         keyboard
%     end

end
end

%% Reassign MLD to LAB if "lith" velocity too slow
% slow_avvlith = 4.56;
% slow_fastlith = 0.05;
% countslow = 0;
labswap = false(Nx,Ny);
% for ix = 1:Nx
% for iy = 1:Ny
%     % don't bother if too far from a seismic station
%     if pt_dist_km(ix,iy) > mindist2sta, continue; end
%     if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
%     % if lith too slow, assign any mld to lab
%     if avvlithsmth(ix,iy) < slow_avvlith ||...
%        fastlithsmth(ix,iy) < slow_fastlith
%         if any(~isnan(zmld_surf_all(ix,iy,:)))
%         zlab_surf_all(ix,iy,:) = zmld_surf_all(ix,iy,:);
%         zmld_surf_all(ix,iy,:) = nan;
%         countslow = countslow+1;
%         labswap(ix,iy) = true;
%         end
%     end
% end
% end
% fprintf('replaced %.0f lab with mld at slow lith\n',countslow)


%% concatenate LAB and MLD surfaces
znvg_surf_all = cat(4,zlab_surf_all,zmld_surf_all);

%% make PRELIM differential NVG depths
znvg_diff_all_prelim = nan(size(znvg_surf_all));
for invg = 1:2
for im = 1:length(methb)
    znvg_diff_all_prelim(:,:,im,invg) =  znvg_surf_all(:,:,im,invg) - nanmean(znvg_surf_all(:,:,im,invg),'all');
end
end


%% QC on the depth picks
% find lonely or outlier points
fprintf('finding lonely/outliers\n')
for invg = 1:2
for im = 1:length(methb)
    fprintf('>> %s-%s\n',nvgs(invg),methb(im))

    
    zzz = znvg_surf_all(:,:,im,invg); 
    fprintf('Nnnan = %.0f\n',sum(~isnan(zzz),'all'))
%     figure(78),clf
%     subplot(131),geogax(gca)
%     surface(longrid,latgrid,zeros(size(longrid)),znvg_surf_all(:,:,im,invg),'linestyle','none')
%     colorbar
    zzz = kill_lonely_outliers(longrid,latgrid,zzz,80,0.7,inf);    
    fprintf('Nnnan = %.0f\n',sum(~isnan(zzz),'all'))
%     subplot(132),geogax(gca)
%     surface(longrid,latgrid,zeros(size(longrid)),a,'linestyle','none')
%     colorbar
    zzz = kill_lonely_outliers(longrid,latgrid,zzz,80,0.9,2);
    fprintf('Nnnan = %.0f\n',sum(~isnan(zzz),'all'))
    zzz = kill_lonely_outliers(longrid,latgrid,zzz,80,0.7,2);
    fprintf('Nnnan = %.0f\n',sum(~isnan(zzz),'all'))
    zzz = kill_lonely_outliers(longrid,latgrid,zzz,80,0.9,1.7);
    fprintf('Nnnan = %.0f\n',sum(~isnan(zzz),'all'))
    
%     subplot(133),geogax(gca)
%     surface(longrid,latgrid,zeros(size(longrid)),b,'linestyle','none')
%     colorbar
%     return
    znvg_surf_all(:,:,im,invg) = zzz;
end
end

% decide on mask - nan out if two method gives nan nvg
% decide on mask - nan out if methods differ in pick by too much
z34diff = [30 10];
z12diff = [50 10];
z_too_muchmean_diff = [80 15];
znvg_nan = false(Nx,Ny,2);
for invg = 1:2  
    count = [0;0;0;0];
for ix = 1:Nx
for iy = 1:Ny
    % too far from a seismic station - nan by default
    if pt_dist_km(ix,iy) > mindist2sta
        znvg_nan(ix,iy,invg) = true;
    % outside of region with good long-period SW data coverage
    elseif ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy))
        znvg_nan(ix,iy,invg) = true;
    % how many methods are nan?
    elseif sum(isnan(znvg_surf_all(ix,iy,:,invg)),3)>=2
        znvg_nan(ix,iy,invg) = true;
        count(1) = count(1)+1;
    
    % methods 3 and 4 differ by too much
    elseif abs(diff(znvg_diff_all_prelim(ix,iy,[4,3],invg)))>z34diff(invg)
    znvg_nan(ix,iy,invg) = true;
         count(2) = count(2)+1;

    % methods 1 and 2 differ by too much
    elseif abs(diff(znvg_diff_all_prelim(ix,iy,[2,1],invg)))>z12diff(invg)
    znvg_nan(ix,iy,invg) = true;
         count(2) = count(2)+1;
   
    % method 3 abs shallower than method 1 abs; shouldn't happen
    elseif diff(znvg_surf_all(ix,iy,[1,3],invg))<0
    znvg_nan(ix,iy,invg) = true;
        count(3) = count(3)+1;
    
    % all methods differ by too much
    elseif sum(abs(demean(squeeze(znvg_diff_all_prelim(ix,iy,:,invg))))) > z_too_muchmean_diff
    znvg_nan(ix,iy,invg) = true;
        count(4) = count(4)+1;
    
    end
end
end
count;
end



% plot ifnan
figure(71);clf
for invg = 1:2
    subplot(1,2,invg)
    pcolor(longrid,latgrid,1+znvg_nan(:,:,invg))
    geogax(gca)
    title("isnan - "+nvgs(invg))
end

%% lab indiv, unsmoothed
figure(912); clf, set(gcf,'pos',[16 477 1497 371])
for im = 1:length(metha)
subplot(1,4,im),hold on
val = znvg_surf_all(:,:,im,1);
val(znvg_nan(:,:,1)) = nan;
surface(longrid,latgrid,zeros(size(longrid)),val,'linestyle','none')
geogax(gca);
title("LAB unsmoothed - "+methb(im))
colormap(brewermap(15,'PuOr')),caxis([70 280])
plot(longrid(labswap),latgrid(labswap),'k.')
if im==4;ax = addcolourbar(gca,gcf); end
end

%% mld indiv, unsmoothed
figure(913); clf, set(gcf,'pos',[16 77 1497 371])
for im = 1:length(metha)
subplot(1,4,im),hold on
val = znvg_surf_all(:,:,im,2);
val(znvg_nan(:,:,2)) = nan;
surface(longrid,latgrid,zeros(size(longrid)),val,'linestyle','none')
geogax(gca);
title("MLD - "+methb(im))
colormap(brewermap(15,'BrBG')),caxis([60 120])
if im==4;ax = addcolourbar(gca,gcf); end
end

%% lab indiv strength
figure(914); clf, set(gcf,'pos',[16 477 1497 371])
for im = 1:length(metha)
subplot(1,4,im),hold on
val = vlab_surf_all(:,:,im);
val(znvg_nan(:,:,1)) = nan;
surface(longrid,latgrid,zeros(size(longrid)),val,'linestyle','none')
geogax(gca);
title("LAB velocity - "+methb(im))
colormap(brewermap(15,'PuOr')),caxis([4.4 4.8]),%caxis([0 5])
plot(longrid(labswap),latgrid(labswap),'k.')
if im==4;ax = addcolourbar(gca,gcf); end
end

%% mld indiv strength
figure(915); clf, set(gcf,'pos',[16 77 1497 371])
for im = 1:length(metha)
subplot(1,4,im),hold on
val = vmld_surf_all(:,:,im);
val(znvg_nan(:,:,2)) = nan;
surface(longrid,latgrid,zeros(size(longrid)),val,'linestyle','none')
geogax(gca);
title("MLD velocity - "+methb(im))
colormap(brewermap(15,'BrBG')),caxis([4.3 4.7]),%caxis([0 5])
if im==4;ax = addcolourbar(gca,gcf); end
end

%% ========================================================
%% Establish mean depths, and some simple stats
znvg_mean = zeros(length(methb),2);
for invg = 1:2 
fprintf('\n%s info\n',upper(nvgs(invg)))
for im = 1:length(methb)
    fprintf('%s method:\n',methb(im))
    avals = znvg_surf_all(:,:,im,invg);
    avals(znvg_nan(:,:,invg)) = nan;

    znvg_mean(im,invg) = mean(avals(~isnan(avals)),'all');
    fprintf('  Mean depth = %.1f km\n',znvg_mean(im,invg))
    fprintf('  Npoints = %.0f\n',sum(~isnan(avals),'all'))
end
end

%% now make confirmed diff maps
znvg_diff_all = nan(size(znvg_surf_all));
for invg = 1:2
for im = 1:length(methb)
    znvg_diff_all(:,:,im,invg) =  znvg_surf_all(:,:,im,invg) - znvg_mean(im,invg);
end
end



%% make map of nvg uncertainty = prop to how much zlab estimates diverge
znvg_std = nan(Nx,Ny,2);
for invg = 1:2    
    znvg_std(:,:,invg) = std(znvg_surf_all(:,:,:,invg),1,3);
end


% feed back - nan out the stds where nan
for invg = 1:2    
    znvg_std_i = znvg_std(:,:,invg);
    znvg_std_i(znvg_nan(:,:,invg)) = nan;
    znvg_std(:,:,invg) = znvg_std_i;
end



% plot uncertainty
figure(76);clf
for invg = 1:2
    subplot(1,2,invg)
    pcolor(longrid,latgrid,znvg_std(:,:,invg))
    geogax(gca)
    title(nvgs(invg))
end

%% See how different estimates vary - abs values
for invg = 1:2
figure(43+invg);clf, set(gcf,'pos', [-200 96 848 770] + 400*invg*[1 0 0 0])
for im1 = 1:3
for im2 = (im1+1):4
subplot(3,3,[im1 + 3*(im2-2)])
hold on
plot([0;300],[0;300],'k--')
s1 = znvg_surf_all(:,:,im1,invg); % 
s2 = znvg_surf_all(:,:,im2,invg);
plot(s1(:),s2(:),'or')
s1(znvg_nan(:,:,invg)) = nan;
s2(znvg_nan(:,:,invg)) = nan;
plot(s1(:),s2(:),'ob')

xlabel(methb(im1));ylabel(methb(im2)); 
if invg == 1
set(gca,'xlim',[50 300],'ylim',[50 300])
elseif invg == 2
set(gca,'xlim',[50 130],'ylim',[50 130])
end

end
end
end
%% See how different estimates vary - diff values
for invg = 1:2
figure(63+invg);clf, set(gcf,'pos', [-200 96 848 770] + 400*invg*[1 0 0 0])
for im1 = 1:3
for im2 = (im1+1):4
subplot(3,3,[im1 + 3*(im2-2)])
hold on
plot([0;300],[0;300],'k--')
s1 = znvg_diff_all(:,:,im1,invg); % 
s2 = znvg_diff_all(:,:,im2,invg);
plot(s1(:),s2(:),'or')
s1(znvg_nan(:,:,invg)) = nan;
s2(znvg_nan(:,:,invg)) = nan;
plot(s1(:),s2(:),'oc')

xlabel(methb(im1));ylabel(methb(im2)); 
if invg == 1
set(gca,'xlim',[-100 100],'ylim',[-100 100])
elseif invg == 2
set(gca,'xlim',[-30 30],'ylim',[-30 30])
end

end
end
end

%% ===============================================
%% do the SMOOTHING
% get smoothed absolute NVG depths
znvg_surfsmooth_all = nan(size(znvg_surf_all));
for invg = 1:2
for im = 1:length(metha)
    for ix = 1:Nx
    for iy = 1:Ny
        if pt_dist_km(ix,iy)>70, continue; end
        if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
        znvg_surfsmooth_all(ix,iy,im,invg) = ...
            gaussian_smoothing_lab(ix,iy,longrid,latgrid,znvg_surf_all(:,:,im,invg),1./znvg_std(:,:,invg).*~znvg_nan(:,:,invg),smthdist);
    end
    end
end
end

% get smoothed differential NVG depths
znvg_diffsmooth_all = nan(size(znvg_diff_all));
for invg = 1:2
for im = 1:length(metha)
    for ix = 1:Nx
    for iy = 1:Ny
        if pt_dist_km(ix,iy)>70, continue; end
        if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
        znvg_diffsmooth_all(ix,iy,im,invg) = ...
            gaussian_smoothing_lab(ix,iy,longrid,latgrid,znvg_diff_all(:,:,im,invg),1./znvg_std(:,:,invg).*~znvg_nan(:,:,invg),smthdist);
    end
    end
end
end

%% Make preferred maps 
% $ cannot simply average, as this assumes that the thickness/breadth of the 
% $ LAB is the same everywhere. If some areas have broader gradients than 
% $ others, then this averaging will mess it all up. This is an argument for 
% $ just using methods "stat" and "v50"
% pin to mean of v50, as a halfway house

% % % use average differential value from all the methods,
% % zlab_pref = znvg_mean(2,1) + ...
% %                 nanmean(znvg_diffsmooth_all(:,:,:,1),3);
% % zmld_pref = znvg_mean(2,2) + ...
% %                 nanmean(znvg_diffsmooth_all(:,:,:,2),3);
% use average differential value from stat and v50 methods,
zlab_pref = znvg_mean(2,1) + ...
                nanmean(znvg_diffsmooth_all(:,:,1:2,1),3);
% zmld_pref = znvg_mean(2,2) + ...
%                 nanmean(znvg_diffsmooth_all(:,:,1:2,2),3);

% preferred depth for MLD is the center of the gradient - this is where RFs
% would show it. The amplitude and Vs, however, is to the v minimum
zmld_pref = znvg_surf_all(:,:,2,2);
vmld_pref = vmld_surf_all(:,:,4);
dvmld_pref = dvmld_surf_all(:,:,4);


% zmld_pref = znvg_surf_all(:,:,4,2);
% vmld_pref = vmld_surf_all(:,:,4);
% dvmld_pref = dvmld_surf_all(:,:,4);


% and now pref but without nan areas
zlab_pref_nonan = zlab_pref; zlab_pref_nonan(znvg_nan(:,:,1)) = nan;

zmld_pref_nonan = zmld_pref; zmld_pref_nonan(znvg_nan(:,:,2)) = nan;


%% compare different estimates to T-based estimate
zlab_T = load('LAB_T1150.mat'); zlab_T = zlab_T.z_lab_Tiso_smth;
figure(66);clf, set(gcf,'pos',[263 813 1271 422])
for im = 1:4
subplot(1,5,im)
hold on
plot([0;300],[0;300],'k--')
s1 = zlab_T;
s2 = znvg_surfsmooth_all(:,:,im,1); % 
plot(s1(:),s2(:),'or')
% s1(znvg_nan(:,:,invg)) = nan;
s2(znvg_nan(:,:,invg)) = nan;
plot(s1(:),s2(:),'o','color',[0.1 0.5 0.6])
xlabel('T-based');ylabel(methb(im)); 
set(gca,'xlim',[70 300],'ylim',[70 300])
end
subplot(1,5,5)
hold on
plot([0;300],[0;300],'k--')
plot(zlab_T(:),zlab_pref_nonan(:),'o','color',[0.1 0.5 0.6])
xlabel('T-based');ylabel('Pref'); 
set(gca,'xlim',[70 300],'ylim',[70 300])


%% Cross sections
figure(5);
lola = [-86.1931   42.6161
        -75.3747   35.2934];
% lola = ginput(2);
[ss,ixs,iys] = points_along_section(lola(1,:),lola(2,:),longrid,latgrid,10);
hold on
plot(lola(:,1),lola(:,2),'pk--','linewidth',3);
%% cross section fig
figure(55);clf,set(gcf,'pos',[62 406 1220 460]); hold on
hexagg = 280;
% dv_ds = 10.;
refV = 4.5;
zlab_T = load('LAB_T1200.mat'); zlab_T = zlab_T.z_lab_Tiso_smth;
for iss = 1:length(ss)
    Vs = squeeze(mgrid3d(ixs(iss),iys(iss),:)) - refV;
    dxv = ss(iss) ;
    ptchVs = Vs;ptchVs(ptchVs<0) = 0;
    % v patch
    patch(([0;ptchVs;0;0]*hexagg)+dxv,depths([1,1:end,end,1]),'r','facealpha',0.2)
    % v profile
    plot((Vs*hexagg)+dxv,depths,'k','linewidth',1.5)
    xline(dxv,':')
    % moho
    plot(dxv + hexagg*[-0.05 0.05],zmoh_surf(ixs(iss),iys(iss))*[1 1],'b','linewidth',2)
    % mld
    plot(dxv + hexagg*[-0.05 0.05],zmld_pref_nonan(ixs(iss),iys(iss))*[1 1],'m','linewidth',2)
%     plot(dxv + (vmld_surf_all(ixs(iss),iys(iss),2)-refV)*hexagg,znvg_surf_all(ixs(iss),iys(iss),2,2),'dm','linewidth',2)
    % lab
    plot(dxv + hexagg*[-0.05 0.05],zlab_pref_nonan(ixs(iss),iys(iss))*[1 1],'r','linewidth',2)
%     plot(dxv + (vlab_surf_all(ixs(iss),iys(iss),1)-refV)*hexagg,znvg_surf_all(ixs(iss),iys(iss),1,1),'dr','linewidth',2)
%     plot(dxv + (vlab_surf_all(ixs(iss),iys(iss),1)-refV)*hexagg,znvg_surf_all(ixs(iss),iys(iss),2,1),'xr','linewidth',2)
%     plot(dxv + (vlab_surf_all(ixs(iss),iys(iss),4)-refV)*hexagg,znvg_surf_all(ixs(iss),iys(iss),4,1),'dr','linewidth',2)
    plot(dxv + hexagg*[-0.05 0.05],zlab_T(ixs(iss),iys(iss))*[1 1],'color',[0.1 0.5 0.6],'linewidth',2)

end
set(gca,'ydir','reverse','ylim',[25 280],'xlim',[0 ss(end)+30])
ylabel('depth (km)','fontsize',20)
xlabel('Distance along section)','fontsize',20)

%% PLOTTING
addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')



%moho
figure(901); clf, hold on
surf(longrid,latgrid,zmoh_surf)
[latbord, lonbord] = borders('states'); % add states map
geogax(gca);
colorbar,colormap(flipud(parula(15))),caxis([30 50])

%% lab pref
figure(902); clf, set(gcf,'pos',[12 302 593 564]);hold on
ZZ = -zlab_pref; ZZ(znvg_nan(:,:,1)) = nan;
surface(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
geogax(gca);
title('LAB - pref')
colorbar,colormap(flipud(brewermap(15,'PuOr'))),caxis(-[260 100])
%% MLD pref
% %mld pref - OLD
% figure(903); clf, hold on
% surface(longrid,latgrid,zeros(size(longrid)),zmld_pref,'linestyle','none')
% geogax(gca)
% title('MLD - pref')
% colorbar,colormap(brewermap(15,'BrBG')),caxis([50 120])

%mld pref - NEW
im = 4; % amplitude of MLD that we care about is vmi method .
zval = zmld_pref_nonan;
% zval(znvg_nan(:,:,2)) = nan;
vval = vmld_pref;
vval(znvg_nan(:,:,2)) = nan;
dvval = dvmld_pref;
dvval(znvg_nan(:,:,2)) = nan;

figure(1903); clf, hold on
scatter(longrid(:),latgrid(:),1000*log10(zval(:)./40).^2,vval(:),'filled','markeredgecolor',0.2*[1 1 1])
geogax(gca);
title('MLD - pref (Vs, km/s)')
colormap(brewermap(15,'BrBG')),caxis([4.35 4.7]),colorbar

figure(1904); clf, hold on
scatter(longrid(:),latgrid(:),1000*log10(zval(:)./40).^2,100*dvval(:),'filled','markeredgecolor',0.2*[1 1 1])
geogax(gca);
title('MLD - pref (dVs, %)')
colormap(flipud(brewermap(15,'PiYG'))),caxis([3 6]),colorbar

figure(1905); clf, hold on
scatter(longrid(:),latgrid(:),1000*log10(zval(:)./40).^2,zval(:),'filled','markeredgecolor',0.2*[1 1 1])
geogax(gca);
title('MLD - pref (depth km)')
colormap(flipud(brewermap(15,'PRGn'))),caxis([70 100]),colorbar


%% diff-NVG maps, smoothed, unsmoothed, all methods
for invg = 1:2
figure(443+invg),clf, set(gcf,'pos',[1 13 1498 853])
if invg == 1, cmap = brewermap(15,'PuOr'); clims = [-80 80]; end
if invg == 2, cmap = brewermap(15,'BrBG'); clims = [-30 30]; end

for im = 1:length(methb)
    subplot(2,4,im);
    ZZ = znvg_diff_all(:,:,im,invg);ZZ(znvg_nan(:,:,invg)) = nan;
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
    geogax(gca);
    colormap(cmap),caxis(clims)
    title("diff-unsmooth-"+methb(im),'fontsize',30);

    subplot(2,4,im+4); hold on
    ZZ = znvg_diffsmooth_all(:,:,im,invg);
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.2)
    ZZ(znvg_nan(:,:,invg)) = nan;
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
    geogax(gca);
    colormap(cmap),caxis(clims)
    title("diff-smooth-",'fontsize',20);
end
end


% look at histograms for range of differential values
figure(446),clf
for im = 1:4
    subplot(2,2,im);
    aa = znvg_diffsmooth_all(:,:,im,1); aa = aa(:);
    hist(aa)
  
    title(metha(im),'fontsize',30);
end


%% saving
if ifsave
LAB_info = struct('longrid',longrid,'latgrid',latgrid,...
                  'zlab_pref',zlab_pref,...
                  'zlab_isnan',znvg_nan(:,:,1),'zlab_wt',1./znvg_std(:,:,1),...
                  'zlab_surf_all',zlab_surf_all,...
                  'zlab_methods',metha,'gauss_smth_dist_km',smthdist,'Vs_smthpts',smthnpts);


MLD_info = struct('longrid',longrid,'latgrid',latgrid,...
                  'zmld_pref',zmld_pref,'vmld_pref',vmld_pref,'dvmld_pref',dvmld_pref,...
                  'zmld_isnan',znvg_nan(:,:,2),'zmld_wt',1./znvg_std(:,:,2),...
                  'zmld_surf_all',zmld_surf_all,...
                  'zmld_methods',metha,'gauss_smth_dist_km',smthdist,'Vs_smthpts',smthnpts);

save('LAB_vgrads','LAB_info','-v7.3')
save('MLD_vgrads','MLD_info','-v7.3')
end

return
%% ========================================================================
%% ========================================================================
%% ========================================================================
%% LAB,MLD nice figure
% three columns, two rows
%  - left shows min NVG (stat) smoothed and differential
%  - right shows max NVG (vmi) smoothed and differential
%  - middle shows pref NVG (mean-v50 + mean(diff_lab)) smoothed and differential
% axpref(row,col)
for invg = 1:2
if invg == 1, cmap = brewermap(15,'PuOr'); end
if invg == 2, cmap = brewermap(15,'BrBG'); end


figure(488 + invg),set(gcf,'position',[512 465 1050 872]),clf
axw = 0.265; axh = 0.4;
axpref(1,1) = axes('position',[0.04 0.52 axw axh]); hold on
axpref(1,2) = axes('position',[0.34 0.52 axw axh]); hold on
axpref(1,3) = axes('position',[0.64 0.52 axw axh]); hold on
axpref(2,1) = axes('position',[0.04 0.07 axw axh]); hold on
axpref(2,2) = axes('position',[0.34 0.07 axw axh]); hold on
axpref(2,3) = axes('position',[0.64 0.07 axw axh]); hold on

% STAT-ABS
ZZ = znvg_surfsmooth_all(:,:,1,invg);
% surf(axpref(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(znvg_surf_all(:,:,1,invg))) = nan;
surf(axpref(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% PREF-ABS
if invg == 1
    ZZ = zlab_pref;
elseif invg == 2
    ZZ = zmld_pref;
end
% surf(axpref(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(znvg_nan(:,:,invg)) = nan;
surf(axpref(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% VMI-ABS
ZZ = znvg_surfsmooth_all(:,:,4,invg);
% surf(axpref(1,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(znvg_surf_all(:,:,4,invg))) = nan;
surf(axpref(1,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% STAT-DIFF
ZZ = znvg_diffsmooth_all(:,:,1,invg);
% surf(axpref(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(znvg_surf_all(:,:,1,invg))) = nan;
surf(axpref(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% PREF-DIFF
if invg == 1
    ZZ = zlab_pref-znvg_mean(2,1);
elseif invg == 2
    ZZ = zmld_pref-znvg_mean(2,2);
end
% surf(axpref(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(znvg_nan(:,:,invg)) = nan;
surf(axpref(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% VMI-DIFF
ZZ = znvg_diffsmooth_all(:,:,4,invg);
% surf(axpref(2,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
ZZ(isnan(znvg_surf_all(:,:,4,invg))) = nan;
surf(axpref(2,3),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')

% adjust colourmaps
if invg == 1, clims1 = [90 250]; clims2 = [-60 60]; end
if invg == 2, clims1 = [70 120]; clims2 = [-20 20]; end
for ia = 1:3
    colormap(axpref(1,ia),cmap);caxis(axpref(1,ia),clims1)
    colormap(axpref(2,ia),cmap);caxis(axpref(2,ia),clims2)
end
% put on geography
for imp = 1:numel(axpref)    
    geogax(axpref(imp));
end
% put on CAMP
a = load('../CAMP_digitize/CAMP_fromGao2000f1_nan.mat');
CAMP = a.a;
for imp = 1:numel(axpref)    
    plot(axpref(imp),CAMP(:,2),CAMP(:,3),'k','Linewidth',1.5);
end

%ticks
set(axpref(:,[2,3]),'yticklabel','')
% titles
titles = {'NVG depth (min)','NVG depth (pref)','NVG depth (max)'};
for ia = 1:3
    title(axpref(1,ia), ...
        regexprep(titles{ia},'NVG',upper(nvgs(invg))), ...
        'fontsize',25,'fontweight','bold')
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

if ifsave
    save2pdf(488 + invg,(upper(nvgs(invg))+'_maps'));
end

end

return

% %% MLD nice figure
% % two columns, two rows
% %  - left shows stat MLD smoothed and differential
% %  - right shows vmi MLD smoothed and differential
% % axpref(row,col)
% 
% figure(490),set(gcf,'position',[1355 465 740 872]),clf
% axw = 0.39; axh = 0.4;
% axprefm(1,1) = axes('position',[0.04 0.52 axw axh]); hold on
% axprefm(1,2) = axes('position',[0.48 0.52 axw axh]); hold on
% axprefm(2,1) = axes('position',[0.04 0.07 axw axh]); hold on
% axprefm(2,2) = axes('position',[0.48 0.07 axw axh]); hold on
% 
% % STAT-ABS
% ZZ = zmld_surfsmooth_all(:,:,1);
% % surf(axprefm(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
% ZZ(isnan(zmld_surf_all(:,:,1))) = nan;
% surf(axprefm(1,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
% 
% % VMI-ABS
% ZZ = zmld_surfsmooth_all(:,:,2);
% % surf(axprefm(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
% ZZ(isnan(zmld_surf_all(:,:,2))) = nan;
% surf(axprefm(1,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
% 
% % STAT-DIFF
% ZZ = zmld_diffsmooth_all(:,:,1);
% % surf(axprefm(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
% ZZ(isnan(zmld_surf_all(:,:,1))) = nan;
% surf(axprefm(2,1),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
% 
% % VMI-DIFF
% ZZ = zmld_diffsmooth_all(:,:,2);
% % surf(axprefm(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.4)
% ZZ(isnan(zmld_surf_all(:,:,2))) = nan;
% surf(axprefm(2,2),longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
% 
% % adjust colourmaps
% for ia = 1:2
%     colormap(axprefm(1,ia),flipud(parula(20))),caxis(axprefm(1,ia),[50 120])
%     colormap(axprefm(2,ia),flipud(turbo(20))),caxis(axprefm(2,ia),[-20 20])
% end
% % put on geobgraphy
% for imp = 1:numel(axprefm)    
%     geogax(axprefm(imp));
% end
% %ticks
% set(axprefm(:,2),'yticklabel','')
% % titles
% titles = {'MLD depth (stat)','MLD depth (vmi)'};
% for ia = 1:2
%     title(axprefm(1,ia),titles{ia},'fontsize',25,'fontweight','bold')
% end
% % colorbars
% cbtitles = {'Absolute depth (km)','Relative depth (km)'};
% for ii = 1:2
%     hcb(ii) = colorbar(axprefm(ii,2),'eastoutside','fontsize',15); drawnow
%     ap = axpos(axprefm(ii,2)); ap(3) = axpos(axprefm(ii,1),3);
%     set(axprefm(ii,2),'position',ap)
%     ylabel(hcb(ii),cbtitles{ii},'fontsize',18,'fontweight','bold')
% end
% set(gcf,'color','w')

if ifsave
    save2pdf(490,'MLD_maps');
end

return

%% Query the map, see the profile
figure(902), pause(0.1)
% figure(5), pause(0.1)
try, delete(hcc); end
drawnow
querypt = ginput(1);
[~,iyq,ixq] = mingrid(abs(longrid-querypt(1)) + abs(latgrid-querypt(2)))
hcc = plot3(longrid(ixq,iyq),latgrid(ixq,iyq),301,'pk','linewidth',2,'markersize',15,'markerfacecolor','r');
Vs = squeeze(mgrid3d(ixq,iyq,:));
[lab_prov,mld_prov] = ...
    LAB_finder(Vs,depths,zmoh_surf(ixq,iyq),smthnpts,1,minmlddv,minlabdv,mldmax)
for jj = 1:3, ha(jj) = subplot(1,3,jj); end
set(ha,'ylim',[40 300]),set(ha(1),'xlim',[4.4 4.8])

%% plot smooth and unsmooth diff labs
for im = 1:4
    subplot(2,2,im); hold on
    ZZ = znvg_diffsmooth_all(:,:,im,1);
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none','FaceAlpha',0.2)
    ZZ(isnan(zlab_surf_all(:,:,im))) = nan;
    surf(longrid,latgrid,zeros(size(longrid)),ZZ,'linestyle','none')
    [latbord, lonbord] = borders('states'); % add states map
    for iplace = 1:length(lonbord)
        line(lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
    end
    xlim([-90 -70]);ylim([27,45]), view(0,90)
    colorbar,colormap(flipud(parula(15))),caxis([-40 40])
    title(metha(im),'fontsize',30);
end

    



%% do a little post-processing
smthnpts = 5;
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


        [lab_prov,mld_prov] = ...
            LAB_finder(Vs1,depths,zmoh_surf(ix,iy),smthnpts,1);
        clone_figure(99,100)
        [lab_prov,mld_prov] = ...
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
zlab_surfsmth_stat(pt_dist_km>mindist2sta | ~geog_inbounds(latgrid,longrid) ) = nan;


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
% function z_smth = gaussian_smoothing_lab(ix,iy,longrid,latgrid,zgrid,wtgrid,smthdist)
% % function to apply a gaussian smoothing function across a surface,
% % smoothing point ix, iy (i.e. grid(ix,iy)) using the points around it.
% % smthdist is the sigma for the gaussian (i.e. 1/4 of the 95% interval
% 
% d2k = 111.1949;
% % find distance IN KM to all other points and use to weight
% dX = abs(longrid - longrid(ix,iy))*d2k*cosd(latgrid(ix,iy)); % account for sphericity in deg2km
% dY = abs(latgrid - latgrid(ix,iy))*d2k;
% dR2 = dX.^2 + dY.^2; % square distance, in km^2
% % ignore all distances greater than 3 sigma
% dR2(dR2 > 9*smthdist.^2) = nan;
% dRwt = exp(-dR2./(2*smthdist.^2));
% 
% nnan = ~isnan(zgrid) & ~isnan(wtgrid) & ~isnan(dR2);
% 
% % calculate the smoothed value
% z_smth = sum(zgrid(nnan).*dRwt(nnan).*wtgrid(nnan))./sum(dRwt(nnan).*wtgrid(nnan));
% end
% 
% function ax = geogax(ax)
%     addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
%     [latbord, lonbord] = borders('states'); % add states map
%     for iplace = 1:length(lonbord)
%         line(ax,lonbord{iplace}',latbord{iplace}','LineWidth',1,'color',0*[1 1 1])
%     end
%     xlim(ax,[-88 -70]);
%     ylim(ax,[25,47]);
%     view(ax,0,90)
%     set(ax,'fontsize',15,'linewidth',2,'box','on','layer','top','Xgrid','off','Ygrid','off')
% end

function ax = addcolourbar(ax,fig)
axsh = get(fig,'Children');
colorbar(ax);
set(ax,'position',[axpos(ax,[1,2]),axpos(axsh(1),[3,4])])
end

function [zmap,ikill] = kill_lonely_outliers(longrid,latgrid,zmap,radius_km,nan_neighbours_thresh_f,nstd_diff_thresh)
% [zmap,ikill] = kill_lonely_outliers(longrid,latgrid,zmap,radius_km,nan_neighbours_thresh_f,nstd_diff_thresh)
%
% function to loop through points and assess how many other non-nan values
% are in their vicinity (test if lonely) and whether they differ
% appreciably from the local average (test if outlier). In each case, the
% test is done relative to all other points within radius_km and the former
% tests whether more than a certain fraction (nonnan_neighbours_f, an
% absolute fraction, not a pct) of the neighbours are nan, while the
% second tests if the point's value differs by more than a certain number
% of local standard deviations (nstd_diff_thresh) from the local average.
% This returns a map with points for which the tests failed nanned-out.
[Nx,Ny] = size(zmap);
ikill = false(Nx,Ny);

d2k = 111.1949;
radius_km2 = radius_km.^2;
aaa = nan(Nx,Ny);
for ix = 1:Nx
for iy = 1:Ny
    if isnan(zmap(ix,iy)),continue; end % no point if already nan!
% find distance IN KM to all other points
dX = abs(longrid - longrid(ix,iy))*d2k*cosd(latgrid(ix,iy)); % account for sphericity in deg2km
dY = abs(latgrid - latgrid(ix,iy))*d2k;
dR2 = dX.^2 + dY.^2; % square distance, in km^2
% find points within radius
ypts = dR2<= radius_km2 & dR2~=0;
% calculate fraction of local nan points
v = zmap(ix,iy);
vv = zmap(ypts);
% if too many neighbours nan, kill point
if sum(isnan(vv))/numel(vv) > nan_neighbours_thresh_f
    ikill(ix,iy) = true;
    continue
end
% assume if we get to this point we have enough non-nan points to produce a
% reasonable std estimate. 
% calculate
vvav = mean(vv(~isnan(vv)));
vvst = std(vv(~isnan(vv)));
if abs(v - vvav)/vvst > nstd_diff_thresh
    ikill(ix,iy) = true;
end
end
end
zmap(ikill) = nan;

end

function [ss,ixs,iys] = points_along_section(lola1,lola2,longrid,latgrid,dss)
% function to make points along a section between lat/lon pairs (1,2).
% Increments of distance in km dss, and then the outputs are ss ([Nx1]),
% the distances along the line and ixs,iys (each [Nx1]) the indices in
% longrid/latgrid of the points closest to the points along the line. nans
% in the ixs/iys if the point closest to the next location is the same as
% the last

[ kmlen, az ] = distance_km(lola1(2),lola1(1),lola2(2),lola2(1));
Ns = ceil(kmlen./dss);
ss = [0:Ns-1]*dss;
ixs = nan(Ns,1);iys = nan(Ns,1);
[lass,loss] = reckon_km(lola1(2),lola1(1),ss,az);
% find closest points to line
for iss = 1:Ns
    [~,iys(iss),ixs(iss)] = mingrid(abs(longrid-loss(iss))+abs(latgrid-lass(iss)));
end

% eliminate repeat points
kill = false(Ns,1);
for iss = 1:Ns
    if any(ixs(iss+1:end)+1i*iys(iss+1:end) == ixs(iss)+1i*iys(iss))
        kill(iss) = true;
    end
end
ixs(kill) = [];
iys(kill) = [];
ss(kill) = [];
Ns = length(ss);

% iss = 1;
% while iss<=length(Ns)
%     find(ixs+1i*iys == ixs(iss)+1i*iys(iss))
%         kill(iss) = true;
%     end
% end
% ixs(kill) = [];
% iys(kill) = [];
% ss(kill) = [];
% Ns = length(ss);

lass = nan(Ns,1);loss = nan(Ns,1);
for iss = 1:Ns
    lass(iss) = latgrid(ixs(iss),iys(iss));
    loss(iss) = longrid(ixs(iss),iys(iss));
end
[ kmlen2, az2 ] = distance_km(lola1(2),lola1(1),lass,loss);
ss = kmlen2.*cosd(az2 - az);
% finally, get rid  of any beyond line ends
kill = false(Ns,1);
kill(ss > kmlen | ss < 0) = true;
ixs(kill) = [];
iys(kill) = [];
ss(kill) = [];


end
















