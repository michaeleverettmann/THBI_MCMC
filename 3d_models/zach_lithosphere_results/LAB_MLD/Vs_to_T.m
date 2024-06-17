% clear

addpath('~/Dropbox/MATLAB/lib/brewermap/')

mindepth = 50; % in km. Ignore depths shallower than this
rhoc = 2750;
rhom = 3300;
g = 9.81;
mindist2sta = 70;
ifsave = true;


%% load in 3D grid
load('surface_colated_b1_V7.mat');
% establish grid size etc.
[Nx,Ny] = size(xgrid);
Nz = length(depths);

load("_distance_pt_to_sta.mat");


%% load in lookup table
Vs_lookup = load('../Vs_to_T/VsTP_lookup_takei_nomelt.mat');


% set up results structures
Tgrid3d = nan(Nx,Ny,Nz);             
etagrid3d = nan(Nx,Ny,Nz);             

for ix = 1:Nx
for iy = 1:Ny
    % don't bother if too far from a seismic station
    if pt_dist_km(ix,iy) > mindist2sta, continue; end
    if ~geog_inbounds(latgrid(ix,iy),longrid(ix,iy)), continue; end    
    % get Vs profile here.
    Vs = squeeze(mgrid3d(ix,iy,:));
    Zmoh = zmoh_surf(ix,iy);
    % loop over depth
    for iz = 1:Nz
        if depths(iz)<mindepth, continue; end
        if depths(iz)<Zmoh, continue; end
        % estimate pressure for this depth
        P_Pa = rhoc*g*Zmoh*1e3 + rhom*g*(depths(iz) - Zmoh)*1e3;
        P_GPa = P_Pa/1e9;
        % set up for linear interpolation
        ipu = find(Vs_lookup.P >= P_GPa,1,'first'); %Vs_lookup.P(ipu)
        ipl = find(Vs_lookup.P < P_GPa,1,'last');  %Vs_lookup.P(ipl)
        f_u = (Vs_lookup.P(ipl) - P_GPa)./diff(Vs_lookup.P([ipu,ipl]));
        f_l = (Vs_lookup.P(ipu) - P_GPa)./diff(Vs_lookup.P([ipl,ipu]));
        % grab Vs(T) vector
        Vs_of_T = f_u*Vs_lookup.Vs_ane(:,ipu) + f_l*Vs_lookup.Vs_ane(:,ipl);
        % calculate eta
        eta_of_T = f_u*Vs_lookup.eta_log10(:,ipu) + f_l*Vs_lookup.eta_log10(:,ipl);
        % interpolate to get T for the observed Vs
        Tgrid3d(ix,iy,iz)   = interp1(Vs_of_T,Vs_lookup.T,Vs(iz)*1e3,'linear','extrap');
        etagrid3d(ix,iy,iz) = interp1(Vs_of_T,eta_of_T,Vs(iz)*1e3,'linear','extrap');

%         % while testing melt effect
%         Vs_of_T_nom = f_u*Vs_lookup2.Vs_ane(:,ipu) + f_l*Vs_lookup2.Vs_ane(:,ipl);
%         T_nom = interp1(Vs_of_T_nom,Vs_lookup.T,Vs(iz)*1e3,'linear','extrap');
%         if abs(T_nom-Tgrid3d(ix,iy,iz)) > 20
%             keyboard
%         end
        
        if isnan(Tgrid3d(ix,iy,iz))
            keyboard
        end
    end
end
end

%% plot some horizontal slices, in temperature
zplts = [80:40:300];
figure(55); clf, set(gcf,'pos',[16 145 1072 832])
for iplt = 1:length(zplts)
zplt = zplts(iplt);
iz = mindex(abs(depths-zplt));

ax(iplt) = subplot(2,3,iplt);
geogax(ax(iplt));
surface(ax(iplt),longrid,latgrid,-ones(Nx,Ny),Tgrid3d(:,:,iz),'linestyle',':')

title(['T (C) at z =',num2str(depths(iz)),' km'])
colormap(flipud(brewermap(15,'RdYlBu')))
caxis([1200 1500])
end
subplot(2,3,6)
% set(gca,'visible','off')
colorbar,colormap(flipud(brewermap(15,'RdYlBu')))
caxis([1200 1600])
set(ax(end),'pos',[axpos(ax(end),[1,2]),axpos(ax(1),[3,4])])
prettyfig(figure(55))

%% plot some horizontal slices, in non-adiabatic T
zplts = [80:40:300];
figure(56); clf, set(gcf,'pos',[16 145 1072 832])
for iplt = 1:length(zplts)
zplt = zplts(iplt);
iz = mindex(abs(depths-zplt));

ax(iplt) = subplot(2,3,iplt);
geogax(ax(iplt));
surface(ax(iplt),longrid,latgrid,-ones(Nx,Ny),Tgrid3d(:,:,iz)-adiabatic_geotherm(depths(iz)),'linestyle',':')

title(['\DeltaT_{ad} (C) at z =',num2str(depths(iz)),' km'])
colormap(flipud(brewermap(15,'RdBu')))
caxis([-150 150])
end
subplot(2,3,6)
% set(gca,'visible','off')
colorbar,colormap(flipud(brewermap(15,'RdBu')))
caxis([-150 150])
set(ax(end),'pos',[axpos(ax(end),[1,2]),axpos(ax(1),[3,4])])
prettyfig(figure(56))

%% plot some horizontal slices, in viscosity
zplts = [80:40:300];
figure(57); clf, set(gcf,'pos',[16 145 1072 832])
for iplt = 1:length(zplts)
zplt = zplts(iplt);
iz = mindex(abs(depths-zplt));

ax2(iplt) = subplot(2,3,iplt);
geogax(ax2(iplt));
surface(ax2(iplt),longrid,latgrid,-ones(Nx,Ny),etagrid3d(:,:,iz),'linestyle',':')

title(['log_{10}\eta (Pa s) at z =',num2str(depths(iz)),' km'])
colormap(gca,(brewermap(15,'Spectral')))
caxis([20 23])
end
subplot(2,3,6)
% set(gca,'visible','off')
colorbar,colormap(gca,(brewermap(15,'Spectral')))
caxis([20 23])
set(ax2(end),'pos',[axpos(ax2(end),[1,2]),axpos(ax2(1),[3,4])])


%% save output
if ifsave
    save('T_map','latgrid','longrid','Tgrid3d','etagrid3d')
    save2jpg(55,'T_depth_map_supp')
    save2jpg(56,'dT_depth_map_supp')

end


%% pretty figures
function figh = prettyfig(figh)
axs = get(figh,'children');
axs = axs(strcmp(get(axs(:),'Type'),'axes'));
for ia = 1:length(axs)
pp = axpos(axs(ia));
% make wider
pp(1) = pp(1)-0.02;
pp(3) = pp(3)+0.04;
% make taller
pp(2) = pp(2)-0.03;
pp(4) = pp(4)+0.06;
set(axs(ia),'pos',pp);
end
set(figh,'color','w')

end


