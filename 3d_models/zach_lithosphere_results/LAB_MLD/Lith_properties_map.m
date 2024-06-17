clear all
%%
% This script does a couple of different integrations of the lithospheric
% velocity profiles to determine the regions where the lithosphere is fast
% and thick. 

%% load in 3D grid
load('surface_colated_b1_V7.mat');
% establish grid size etc.
[Nx,Ny] = size(xgrid);
Nz = length(depths);

smthdist = 35;
mindist2sta = 70;

vfastcutoff = 4.5;
zfastcutoff = 250;
ref_rgn = 'cont_LProt'; %cont_Arch cont_EProt cont_LProt cont_Phan

ifsave = true;

addpath('/Users/zeilon/Dropbox/MATLAB/lib/borders_v3.1.2/borders')
addpath('~/Dropbox/MATLAB/lib/brewermap/')
addpath('../model_process_plot/')

load("_distance_pt_to_sta.mat");
zmoh_surf(pt_dist_km>mindist2sta) = nan;


% % try plotting 3D contour
% LO3 = repmat(longrid,[1,1,Nz]);
% LA3 = repmat(latgrid,[1,1,Nz]);
% ZZ3 = repmat(reshape(depths,[1,1,Nz]),[Nx,Ny,1]);
% F46 = isosurface(LO3,LA3,ZZ3,mgrid3d,4.6);
% F47 = isosurface(LO3,LA3,ZZ3,mgrid3d,4.7);
% figure(434),clf, hold on
% p = patch(F46);
% p2 = patch(F47);
% p.FaceColor = 'red';  p2.FaceColor = 'blue';
% p.EdgeColor = 'none'; p2.EdgeColor = 'none';
% p.FaceAlpha = 0.2;    p2.FaceAlpha = 0.2;
% set(gca,'zdir','reverse')
% % colormap(prism(28))
% % camup([1 0 0 ]); campos([25 -55 5]) 


%% Map of how much of the lithosphere is "fast"
fastlith = nan(Nx,Ny);

for ix = 1:Nx
for iy = 1:Ny
    if pt_dist_km(ix,iy)>70, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    vv = squeeze(mgrid3d(ix,iy,:));
    vv = vv-vfastcutoff;
    vv(vv<0)=0;
    vv(depths>zfastcutoff) = 0;
    vv(depths<zmoh_surf(ix,iy)) = 0;
    fastlith(ix,iy) = trapz(depths,vv)./(zfastcutoff-zmoh_surf(ix,iy));
end
end

[fastlith] = kill_lonely_outliers(longrid,latgrid,fastlith,150,0.7,3);

% now a smooth version
fastlithsmth = nan(Nx,Ny);
for ix = 1:Nx
for iy = 1:Ny
    if pt_dist_km(ix,iy)>70, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    fastlithsmth(ix,iy) = ...
        gaussian_smoothing_lab(ix,iy,longrid,latgrid,fastlith,ones(Nx,Ny),smthdist);
end
end
fastlith(pt_dist_km>mindist2sta) = nan;
fastlithsmth(pt_dist_km>mindist2sta) = nan;

% create reference value - integrate SEMUCB cont_LProt profile
nasdatapath = '/Volumes/data/';
addpath([nasdatapath,'models_seismic/SEMum2_avg_VS/']);
a = SEMum2_avgprofiles(0,[nasdatapath,'models_seismic/SEMum2_avg_VS/']);
rgns = fieldnames(a); rgns = setdiff(rgns,'Z');
fastlith_rgn = zeros(length(rgns),1);
for irgn = 1:length(rgns)
    vv = a.(rgns{irgn});
    vv = vv-vfastcutoff;
    vv(vv<0)=0;
    vv(a.Z>zfastcutoff) = 0;
    ivv = (a.Z >= 50 & a.Z <=zfastcutoff & ~isnan(vv));
    fastlith_rgn(irgn) = trapz(a.Z(ivv),vv(ivv))./(max(a.Z(ivv)) - min(a.Z(ivv)));
end

fastlith_ref = fastlith_rgn(strcmp(rgns,ref_rgn));

% normalized fastlith
fastlith_norm = fastlith./fastlith_ref;
fastlithsmth_norm = fastlithsmth./fastlith_ref;

% plot it
figure(5),clf,hold on
surf(longrid,latgrid,-ones(Nx,Ny),fastlithsmth_norm,'facealpha',0.9)
% put a contour on
contour(longrid,latgrid,fastlithsmth_norm,1*[1 1],'linewidth',1.5,'color','k')
contour(longrid,latgrid,fastlithsmth_norm,1.5*[1 1],'linewidth',2,'color','k')
contour(longrid,latgrid,fastlithsmth_norm,2*[1 1],'linewidth',2.5,'color','k')

geogax(gca);
colorbar,colormap(flipud(turbo(15)))
title(['FastLith (relative to ',regexprep(ref_rgn,'_','-'),')'])

%% Map of average lithospheric velocity
avvlith = nan(Nx,Ny);
for ix = 1:Nx
for iy = 1:Ny
    if pt_dist_km(ix,iy)>70, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    vv = squeeze(mgrid3d(ix,iy,:));
    iz = depths<=zfastcutoff & depths>=zmoh_surf(ix,iy)+5;
    avvlith(ix,iy) = trapz(depths(iz),vv(iz))./(max(depths(iz)) - min(depths(iz)));
end
end

[avvlith] = kill_lonely_outliers(longrid,latgrid,avvlith,150,0.6,2);

% now a smooth version
avvlithsmth = nan(Nx,Ny);
for ix = 1:Nx
for iy = 1:Ny
    if pt_dist_km(ix,iy)>70, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    avvlithsmth(ix,iy) = ...
        gaussian_smoothing_lab(ix,iy,longrid,latgrid,avvlith,ones(Nx,Ny),smthdist);
end
end
avvlith(pt_dist_km>mindist2sta) = nan;
avvlithsmth(pt_dist_km>mindist2sta) = nan;

% plot it
figure(6),clf,hold on
surf(longrid,latgrid,-ones(Nx,Ny),avvlithsmth,'facealpha',0.9)
% put a contour on
contour(longrid,latgrid,avvlithsmth,4.52*[1 1],'linewidth',1.5,'color','k')
contour(longrid,latgrid,avvlithsmth,4.57*[1 1],'linewidth',2,'color','k')
contour(longrid,latgrid,avvlithsmth,4.62*[1 1],'linewidth',2.5,'color','k')

geogax(gca);
colorbar,colormap(flipud(parula(15))), caxis([4.48 4.65])
title('Mean UM Vs (km/s)')

%% SAVE
if ifsave
    save('fastlith_map','longrid','latgrid','fastlithsmth_norm','avvlithsmth','ref_rgn')
end

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
















