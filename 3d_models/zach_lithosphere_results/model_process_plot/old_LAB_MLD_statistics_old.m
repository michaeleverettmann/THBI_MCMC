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


% fit surfaces
Ffastlithsmth_norm = scatteredInterpolant(longrid(:),latgrid(:),fastlithsmth_norm(:));
Favvlithsmth = scatteredInterpolant(longrid(:),latgrid(:),avvlithsmth(:));
Fz_lab_Tiso_smth = scatteredInterpolant(longrid(:),latgrid(:),z_lab_Tiso_smth(:));
Fzmld_pref = scatteredInterpolant(longrid(:),latgrid(:),MLD_info.zmld_pref(:));
Fvmld_pref = scatteredInterpolant(longrid(:),latgrid(:),MLD_info.vmld_pref(:));
Fdvmld_pref = scatteredInterpolant(longrid(:),latgrid(:),MLD_info.dvmld_pref(:));

% prep output velocity profiles
vz_rgn = zeros(length(z_med),4);
T_at150 = zeros(1,4);

fprintf('=============================================\n')
%% cycle through regions
regions_name = ["Margin","Craton","Grenville","NE Anomaly"];
for ir = [2,3,1,4]
    % find stations in this region
    ista_rgn = find(stas_info.region==ir);
    Nsta = length(ista_rgn);
    % set up results for averages
    zlab_rgn = nan(Nsta,1);
    fastlith_rgn = nan(Nsta,1);
    avv_rgn = nan(Nsta,1);
    zmld_rgn = nan(Nsta,1);
    vmld_rgn = nan(Nsta,1);
    dvmld_rgn = nan(Nsta,1);
    % loop thorough stas
    for is = 1:Nsta
        % accumulate Vs(z) profile
        vz_rgn(:,ir) = vz_rgn(:,ir) + v_med(ista_rgn(is),:)';
        % find station location, 
        slat = stas_info.lat(ista_rgn(is));
        slon = stas_info.lon(ista_rgn(is));
        [d,iys,ixs] = mingrid(abs(longrid-slon)+abs(latgrid-slat));
        % query fit surfaces
        zlab_rgn(is) = Fz_lab_Tiso_smth(slon,slat);
        fastlith_rgn(is) = Ffastlithsmth_norm(slon,slat);
        avv_rgn(is) = Favvlithsmth(slon,slat);
        if d<0.8 && ~MLD_info.zmld_isnan(ixs,iys)
            zmld_rgn(is) = Fzmld_pref(slon,slat);
            vmld_rgn(is) = Fvmld_pref(slon,slat);
            dvmld_rgn(is) = Fdvmld_pref(slon,slat);
        end
        % grab temperature
        T_at150(ir) = T_at150(ir) + Tgrid3d(ixs,iys,a.depths==150);
    end
    vz_rgn(:,ir) = vz_rgn(:,ir)./Nsta;
    T_at150(ir) = T_at150(ir)./Nsta;
    figure(22),clf
    hist(dvmld_rgn)
    fprintf('\n%s region:\n',regions_name(ir))
    fprintf('Z-LAB    = %4.0f %s %4.0f km \t (N = %2.0f)\n',nanmean(zlab_rgn),char(177),nanstd(zlab_rgn),sum(~isnan(zlab_rgn)))
    fprintf('fastlith = %4.2f %s %4.2f  \t (N = %2.0f)\n',nanmean(fastlith_rgn),char(177),nanstd(fastlith_rgn),sum(~isnan(fastlith_rgn)))
    fprintf('avvlith  = %4.2f %s %4.2f km/s  \t (N = %2.0f)\n',nanmean(avv_rgn),char(177),nanstd(avv_rgn),sum(~isnan(avv_rgn)))
    fprintf('Z-MLD    = %4.0f %s %4.0f km \t (N = %2.0f)\n',nanmean(zmld_rgn),char(177),nanstd(zmld_rgn),sum(~isnan(zmld_rgn)))
    fprintf('V-MLD    = %4.2f %s %4.2f km/s \t (N = %2.0f)\n',nanmean(vmld_rgn),char(177),nanstd(vmld_rgn),sum(~isnan(vmld_rgn)))
    fprintf('dV-MLD   = %4.2f %s %4.2f %% \t (N = %2.0f)\n',100*nanmean(dvmld_rgn),char(177),100*nanstd(dvmld_rgn),sum(~isnan(dvmld_rgn)))
    fprintf('T(150km) = %4.0f C\n', T_at150(ir))

end

figure(23);clf
plot(vz_rgn,z_med,'linewidth',2); 
set(gca,'ydir','reverse','xlim',[4.4 4.8],'ylim',[50 280])
legend(regions_name)

%% query lith thickness drop across appalachian front
pt1 = [-81,38];% se WV, just W of HA
pt2 = [-79.5,37.5]; % just S of HA
dz = Fz_lab_Tiso_smth(pt1(1),pt1(2)) - Fz_lab_Tiso_smth(pt2(1),pt2(2));
dh = distance_km(pt1(2),pt1(1),pt2(2),pt2(1));
fprintf('Difference in thickness of Lith across AF:\n')
fprintf('%.0f km\n',dz)
fprintf('Over horizontal distance of %.0f km\n',dh)
fprintf('(dip of %.0f)\n',atand(dz/dh))

    

