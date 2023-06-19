clear all
%% load in 3D grid
load('../surface_colated_b1_V7.mat');
% establish grid size etc.
[Nx,Ny] = size(xgrid);
Nz = length(depths);

% %% test by plotting moho
figure(901); clf
surface(xgrid,ygrid,zmoh_surf)
colorbar,colormap(flipud(parula))

%% loop over all the vertical profiles and find lab/mld
zlab1_surf = nan(size(zmoh_surf));
zlab2_surf = nan(size(zmoh_surf));
zmld1_surf = nan(size(zmoh_surf));
zmld2_surf = nan(size(zmoh_surf));

for ix = 1:Nx
for iy = 1:Ny
    Vs = squeeze(mgrid3d(ix,iy,:));
    try
        [lab_provisional,mld_provisional] = ...
            LAB_finder(Vs,depths,zmoh_surf(ix,iy),10,0);

        if ~isempty(lab_provisional.z_stat) &&~isnan(lab_provisional.z_stat)
            zlab1_surf(ix,iy) = lab_provisional.z_stat;
        end
        if ~isempty(lab_provisional.z_vmi) && ~isnan(lab_provisional.z_vmi)
            zlab2_surf(ix,iy) = lab_provisional.z_vmi;
        end
        if ~isempty(mld_provisional.z_stat) && ~isnan(mld_provisional.z_stat)
            zmld1_surf(ix,iy) = mld_provisional.z_stat;  
        end
          if ~isempty(mld_provisional.z_vmi) && ~isnan(mld_provisional.z_vmi)
            zmld2_surf(ix,iy) = mld_provisional.z_vmi;  
        end
    catch
        fprintf('Some failure for [%.2f N, %.2f E]\n',latgrid(ix,iy),longrid(ix,iy))
    end
end
end

figure(902); clf
surface(xgrid,ygrid,zlab2_surf)
colorbar,colormap(flipud(turbo))

figure(903); clf
surface(xgrid,ygrid,zmld2_surf)
colorbar,colormap(flipud(turbo)),caxis([50 120])

