%% print some statistics for the LAB and MLD values within the regions
clear all
% load things
load('stas_info.mat');
load('stas_collated.mat','v_med','z_med');
load('../LAB_MLD/LAB_T1150.mat')
load('../LAB_MLD/MLD_vgrads.mat')
load('../LAB_MLD/fastlith_map.mat')
a = load('../LAB_MLD/surface_colated_b1_V7.mat');
load('../LAB_MLD/T_map.mat')

ifsave = true;

%% cycle through lat/lon points, determining region
[Nx,Ny] = size(latgrid);
rgngrid = nan(Nx,Ny);
for ix = 1:Nx
for iy = 1:Ny
    if isnan(z_lab_Tiso_smth(ix,iy)),continue; end
    d = distance_km(latgrid(ix,iy),longrid(ix,iy),stas_info.lat,stas_info.lon);
    okd = d<100;
    modergn = mode(stas_info.region(okd));
    if length(modergn)==1
        rgngrid(ix,iy) = modergn;
    else
        keyboard
    end
end
end
% map it
figure(21),clf, hold on
geogax(gca);
scatter(longrid(:),latgrid(:),40,rgngrid(:)+1,'filled')




% % fit surfaces
% Ffastlithsmth_norm = scatteredInterpolant(longrid(:),latgrid(:),fastlithsmth_norm(:));
% Favvlithsmth = scatteredInterpolant(longrid(:),latgrid(:),avvlithsmth(:));
Fz_lab_Tiso_smth = scatteredInterpolant(longrid(:),latgrid(:),z_lab_Tiso_smth(:));
% Fzmld_pref = scatteredInterpolant(longrid(:),latgrid(:),MLD_info.zmld_pref(:));
% Fvmld_pref = scatteredInterpolant(longrid(:),latgrid(:),MLD_info.vmld_pref(:));
% Fdvmld_pref = scatteredInterpolant(longrid(:),latgrid(:),MLD_info.dvmld_pref(:));

% prep output velocity profiles
vz_rgn = zeros(length(a.depths),4);
T_at150 = zeros(1,4);

fig_L = figure(22);clf, set(fig_L,'Position',[75 288 1130 713],'color','w');
fig_M = figure(23);clf, set(fig_M,'Position',[141 109 1130 540],'color','w')


fprintf('=============================================\n')
%% cycle through regions
regions_name = ["Margin","Craton","Grenville","NE Anomaly"];
irplt =0;
for ir = [2,3,1,4]
    irplt = irplt+1;
    % consider points in this region
    % set up results for averages
    ok_L = rgngrid==ir;
    ok_M = ok_L & ~MLD_info.zmld_isnan;
    zlab_rgn(ir) = median(z_lab_Tiso_smth(ok_L));
    zlab_rgn_std(ir) = nanstd(z_lab_Tiso_smth(ok_L));
    fastlith_rgn(ir) = nanmedian(fastlithsmth_norm(ok_L));
    fastlith_rgn_std(ir) = nanstd(fastlithsmth_norm(ok_L));
    avv_rgn(ir) = nanmedian(avvlithsmth(ok_L));
    avv_rgn_std(ir) = nanstd(avvlithsmth(ok_L));
    zmld_rgn(ir) = nanmedian(MLD_info.zmld_pref(ok_M));
    zmld_rgn_std(ir) = nanstd(MLD_info.zmld_pref(ok_M));
    vmld_rgn(ir) = nanmedian(MLD_info.vmld_pref(ok_M));
    vmld_rgn_std(ir) = nanstd(MLD_info.vmld_pref(ok_M));
    dvmld_rgn(ir) = nanmedian(MLD_info.dvmld_pref(ok_M));
    dvmld_rgn_std(ir) = nanstd(MLD_info.dvmld_pref(ok_M));
    % grab temperature
    Tgrid3d_z = Tgrid3d(:,:,a.depths==150);
    T_at150(ir) = median(Tgrid3d_z(ok_L));
    
    for id = 1:length(a.depths)
        mslc = a.mgrid3d(:,:,id);
        vz_rgn(id,ir) = median(mslc(ok_L));
    end
    vz_rgn(:,ir) = smooth(vz_rgn(:,ir),MLD_info.Vs_smthpts);

    % plot histograms
    % MLD values
    figure(fig_M)
    
    subplotij(3,4,1,irplt); cla, hold on
    histogram(MLD_info.zmld_pref(ok_M),[60:4:110],'normalization','pdf')
    xline(zmld_rgn(ir),'--r','linewidth',1.5)
    title("Region "+regions_name(ir))
    xlabel('MLD depth (km)'),xlim([60 110]),set(gca,'box','on','linewidth',1.5)
    
    subplotij(3,4,2,irplt); cla, hold on
    histogram(MLD_info.vmld_pref(ok_M),[4.3:0.02:4.8],'normalization','pdf')
    xline(vmld_rgn(ir),'--r','linewidth',1.5)
    xlabel('MLD Vs (km/s)'),xlim([4.3 4.8]),set(gca,'box','on','linewidth',1.5)
    
    subplotij(3,4,3,irplt); cla, hold on
    histogram(100*MLD_info.dvmld_pref(ok_M),[0.5:0.2:8],'normalization','pdf')
    xline(100*dvmld_rgn(ir),'--r','linewidth',1.5)
    xlabel('MLD dVs (%)'),xlim([0 8]),set(gca,'box','on','linewidth',1.5)

    % lithosphere values
    figure(fig_L)
    
    subplotij(4,4,1,irplt); cla, hold on
    histogram(z_lab_Tiso_smth(ok_L),[80:10:250],'normalization','pdf')
    xline(zlab_rgn(ir),'--r','linewidth',1.5)
    title("Region "+regions_name(ir))
    xlabel('LAB depth (km)'),xlim([80 260]),set(gca,'box','on','linewidth',1.5)
    
    subplotij(4,4,2,irplt); cla, hold on
    histogram(fastlithsmth_norm(ok_L),[0:0.1:2.1],'normalization','pdf')
    xline(fastlith_rgn(ir),'--r','linewidth',1.5),set(gca,'box','on','linewidth',1.5)
    xlabel('fastlith'),xlim([0.0 2.2])
    
    subplotij(4,4,3,irplt); cla, hold on
    histogram(avvlithsmth(ok_L),[4.3:0.02:4.8],'normalization','pdf')
    xline(avv_rgn(ir),'--r','linewidth',1.5)
    xlabel('avV-lith (km/s)'),xlim([4.45 4.7]),set(gca,'box','on','linewidth',1.5)
    
    subplotij(4,4,4,irplt); cla, hold on
    histogram(Tgrid3d_z(ok_L),[600:50:1250],'normalization','pdf')
    xline(T_at150(ir),'--r','linewidth',1.5)
    xlabel('T at 150 km (C)'),xlim([600 1300]),set(gca,'box','on','linewidth',1.5)


    % print out key values
    fprintf('\n%s region:\n',regions_name(ir))
    fprintf('Z-Moh    = %4.0f km \t (N = %2.0f)\n',nanmedian(a.zmoh_surf(ok_L)),sum(ok_L,'all'))
    fprintf('Z-LAB    = %4.0f %s %4.0f km \t (N = %2.0f)\n',zlab_rgn(ir),char(177),zlab_rgn_std(ir),sum(ok_L,'all'))
    fprintf('fastlith = %4.2f %s %4.2f  \t (N = %2.0f)\n',fastlith_rgn(ir),char(177),fastlith_rgn_std(ir),sum(ok_L,'all'))
    fprintf('avvlith  = %4.2f %s %4.2f km/s  \t (N = %2.0f)\n',avv_rgn(ir),char(177),avv_rgn_std(ir),sum(ok_L,'all'))
    fprintf('Z-MLD    = %4.0f %s %4.0f km \t (N = %2.0f)\n',zmld_rgn(ir),char(177),zmld_rgn_std(ir),sum(ok_M,'all'))
    fprintf('V-MLD    = %4.2f %s %4.2f km/s \t (N = %2.0f)\n',vmld_rgn(ir),char(177),vmld_rgn_std(ir),sum(ok_M,'all'))
    fprintf('dV-MLD   = %4.2f %s %4.2f %% \t (N = %2.0f)\n',100*dvmld_rgn(ir),char(177),100*dvmld_rgn_std(ir),sum(ok_M,'all'))
    fprintf('T(150km) = %4.0f C\n', T_at150(ir))

end

figure(25);clf
plot(vz_rgn,a.depths,'linewidth',2); 
set(gca,'ydir','reverse','xlim',[4.4 4.8],'ylim',[50 280])
legend(regions_name)

% save
if ifsave
    save2jpg(fig_L,'Region_stats_Lith');
    save2jpg(fig_M,'Regions_stats_MLD');
end

%% query lith thickness drop across appalachian front
pt1 = [-81,38];% se WV, just W of HA
pt2 = [-79.5,37.5]; % just S of HA
dz = Fz_lab_Tiso_smth(pt1(1),pt1(2)) - Fz_lab_Tiso_smth(pt2(1),pt2(2));
dh = distance_km(pt1(2),pt1(1),pt2(2),pt2(1));
fprintf('Difference in thickness of Lith across AF:\n')
fprintf('%.0f km\n',dz)
fprintf('Over horizontal distance of %.0f km\n',dh)
fprintf('(dip of %.0f)\n',atand(dz/dh))

%% MLD values
fprintf('\nMLD info  \n')
% SAA
SAA_lola = [-84.9,35];
d2A = 50;
iaa = distance_km(SAA_lola(2),SAA_lola(1),latgrid(:),longrid(:)) < d2A;
fprintf('SAA-MLD:  ')
fprintf('Depth = %3.0f km, V = %4.2f km/s, dV = %3.1f %% \n',...
    median(MLD_info.zmld_pref(iaa)),median(MLD_info.vmld_pref(iaa)),100*median(MLD_info.dvmld_pref(iaa)))
% CAA
CAA_lola = [-79.95,38.08];
d2A = 50;
iaa = distance_km(CAA_lola(2),CAA_lola(1),latgrid(:),longrid(:)) < d2A;
fprintf('CAA-MLD:  ')
fprintf('Depth = %3.0f km, V = %4.2f km/s, dV = %3.1f %% \n',...
    median(MLD_info.zmld_pref(iaa)),median(MLD_info.vmld_pref(iaa)),100*median(MLD_info.dvmld_pref(iaa)))
% NAA
NAA_lola = [-72.4,43.4];
d2A = 100;
iaa = distance_km(NAA_lola(2),NAA_lola(1),latgrid(:),longrid(:)) < d2A;
fprintf('NAA-MLD:  ')
fprintf('Depth = %3.0f km, V = %4.2f km/s, dV = %3.1f %% \n',...
    nanmedian(MLD_info.zmld_pref(iaa)),nanmedian(MLD_info.vmld_pref(iaa)),100*nanmedian(MLD_info.dvmld_pref(iaa)))



%% crustal thickness vs. lithospheric thickness at AF
addpath('../model_process_plot')
AF = appalachian_front;
% loop through all points in model, find closest distance to AF
d2AF_km = zeros(size(latgrid));
AFclosest_lola = zeros(Nx,Ny,2);
for ix = 1:Nx
for iy = 1:Ny
    d2AF_i = distance_km(latgrid(ix,iy),longrid(ix,iy),AF(:,2),AF(:,1));
    [d2AF_km(ix,iy),ii] = min(d2AF_i);
    AFclosest_lola(ix,iy,1:2) = AF(ii,:);
end
end
mind2AF = 140;
lonlims = [-84.5 -75.5];
d2AF_km(d2AF_km>mind2AF) = nan;
d2AF_km(longrid < min(lonlims)) = nan;
d2AF_km(longrid > max(lonlims)) = nan;
%figure;pcolor(d2AF_km')

% now find average crustal thickness in this zone
figure(89),clf, set(gcf,'Position',[45 549 1318 317])
subplot(131)
scatter(longrid(~isnan(d2AF_km)),latgrid(~isnan(d2AF_km)),50,a.zmoh_surf(~isnan(d2AF_km)),'filled')
title('Z=moh'),colorbar, colormap(gca,parula)
subplot(132)
scatter(longrid(~isnan(d2AF_km)),latgrid(~isnan(d2AF_km)),50,z_lab_Tiso_smth(~isnan(d2AF_km)),'filled')
title('Z-LAB (T1150)'),colorbar, colormap(gca,flipud(parula))
subplot(133)
plot(a.zmoh_surf(~isnan(d2AF_km)),z_lab_Tiso_smth(~isnan(d2AF_km)),'o')
xlabel('Z-Moh'),ylabel('Z-LAB')

fprintf('Av Moho depth: %4.1f km\n',median(a.zmoh_surf(~isnan(d2AF_km))))
fprintf('Av LAB depth:  %4.0f km\n',median(z_lab_Tiso_smth(~isnan(d2AF_km))))



    

