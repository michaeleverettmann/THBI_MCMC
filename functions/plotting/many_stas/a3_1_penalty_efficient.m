function [penalty] = a3_1_penalty_efficient(vgrid,...
    pdf_terp, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay,...
    nearesti, weighti, min_mm_terp, dmm_di, nmm, nsta); 

dvdx2 = ((vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2); 
dvdy2 = ((vgrid(:,1:end-2) - 2*vgrid(:,2:end-1) + vgrid(:,3:end))./dy2); 
roughness = sum(dvdx2 .^ 2, 'all') + sum(dvdy2 .^ 2, 'all'); % Using very simple roughness for now. For speed Should change to the below line. 2022.09.27. 
roughness = roughness * rough_scale; 

% vsta = interpn(xgrid, ygrid, vgrid, stax, stay, 'linear'); % Modelled velocity at each station 

% Interpolate v at stations. All the spatial work is already done. 
% nearesti and weightsi come from "p_prep_grid_to_sta_interp"
vsta = sum(vgrid(nearesti) .* weighti,2); 

% Interpolate v at each station to a pdf. 
pv_mod = p_mm_to_pdf_dmdi(vsta, pdf_terp, min_mm_terp, dmm_di, nmm, nsta); 

penalty = - sum(pv_mod); % Lower penalty is higher probability. 
penalty = penalty + roughness; 

end