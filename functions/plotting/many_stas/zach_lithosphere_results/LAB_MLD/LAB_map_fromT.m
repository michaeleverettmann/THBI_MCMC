% clear

addpath('~/Dropbox/MATLAB/lib/brewermap/')
addpath('../model_process_plot/')

mindepth = 50; % in km. Ignore depths shallower than this
mindist2sta = 70;
T_lab = 1150;
intTval_norm_thresh = 0.5;
z_minlab = 70;
z_maxlab = 280;
smthdist = 45;
ifsave = false;


%% load in 3D grid
load('surface_colated_b1_V7.mat');
% establish grid size etc.
[Nx,Ny] = size(xgrid);
Nz = length(depths);

load("_distance_pt_to_sta.mat");


%% load in T maps
load('T_map.mat');


%% now depth to first passing of isotherm
% set up results structures
ifplot = false;
z_lab_Tiso = nan(Nx,Ny);             
for ix = 1:Nx
for iy = 1:Ny
    % don't bother if too far from a seismic station
    if pt_dist_km(ix,iy) > mindist2sta, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    % find isotherm crossing depth(s)
    Tz = squeeze(Tgrid3d(ix,iy,:));
    [~,z_try] = crossing(Tz,depths,T_lab);
    z_try(z_try<z_minlab) = [];
    z_try(z_try>z_maxlab) = [];
    if isempty(z_try),continue; end
    intT_to_pick = 1; % default is just shallowest one

    % now need a way to discriminate if crossing isotherm due to MLD (and
    % possible melt or exogenous mineralogy) or if as real LAB. If the
    % former, expect lots of low temperatures again below the z_try
    % candidate point. If the latter, then almost all of the low
    % temperatures are above this point. 
    % So integrate the temperatures beneath T_lab from z_try to surface.
    % What fraction is found at each depth?
    
    if length(z_try)>1
        Tint = T_lab-Tz; % so, positive if colder
        Tint(Tint<0) = 0; % ignore slow temperatures altogether
        intTval = zeros(size(z_try));
        for itry = 1:length(z_try)
            d2int = (depths(:) >= mindepth) & (depths(:) <= z_try(itry)) & ~isnan(Tint(:));
            intTval(itry) = trapz(depths(d2int),Tint(d2int));
        end
        intTval_norm = intTval./max(intTval);
        intT_to_pick = find(intTval_norm>intTval_norm_thresh,1,"first");

        if ifplot
            figure(564); clf, hold on
            plot(Tz,depths,'linewidth',1.5);
            xline(T_lab,'--r')
            plot(T_lab,z_try,'dr')
            plot(T_lab,z_try(intT_to_pick),'ok','markersize',10,'MarkerFaceColor','k')
            set(gca,'ydir','reverse','box','on','linewidth',1.5),ylabel('Z (km)'),xlabel('T (C)');
            pause
        end
    end
    

%     if length(z_try)>1,keyboard; end
    z_lab_Tiso(ix,iy) = z_try(intT_to_pick);
%     z_lab_Tiso(ix,iy) = median(z_try);
end
end
% now smooth
z_lab_Tiso_smth = nan(Nx,Ny);
for ix = 1:Nx
for iy = 1:Ny
    if pt_dist_km(ix,iy) > mindist2sta, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end
    z_lab_Tiso_smth(ix,iy) = gaussian_smoothing_lab(ix,iy,longrid,latgrid,z_lab_Tiso,1,smthdist);
end
end

%% PLOT
figure(57); clf, set(gcf,'pos',[16 292 589 574])
geogax(gca); hold on
surface(longrid,latgrid,-ones(Nx,Ny),z_lab_Tiso_smth,'linestyle','none')

title(['Z_{LAB} (km) for T = ',num2str(T_lab),' C'])
colormap(brewermap(15,'PuOr'))
caxis([z_minlab 260])
colorbar

% put on CAMP
a = load('../CAMP_digitize/CAMP_fromGao2000f1_nan.mat');
CAMP = a.a;
plot(CAMP(:,2),CAMP(:,3),'k','Linewidth',2);
plot(CAMP(:,2),CAMP(:,3),'g','Linewidth',1);


%% save output
if ifsave
    save(['LAB_T',num2str(T_lab)],'z_lab_Tiso_smth','latgrid','longrid','z_lab_Tiso','smthdist','intTval_norm_thresh')
end




%% subfunctions
