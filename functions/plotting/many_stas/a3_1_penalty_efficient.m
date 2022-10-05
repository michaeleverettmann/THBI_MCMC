function [penalty] = a3_1_penalty_efficient(vgrid,...
    pdfs, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay); 

dvdx2 = ((vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2); 
dvdy2 = ((vgrid(:,1:end-2) - 2*vgrid(:,2:end-1) + vgrid(:,3:end))./dy2); 
roughness = sum(dvdx2 .^ 2, 'all') + sum(dvdy2 .^ 2, 'all'); % Using very simple roughness for now. For speed Should change to the below line. 2022.09.27. 
roughness = roughness * rough_scale; 

vsta = interpn(xgrid, ygrid, vgrid, stax, stay, 'linear'); % Modelled velocity at each station 

% ista = 1; 
pv_mod = zeros(size(stax)); 
for ista=1:length(pv_mod); 
    pv_mod(ista) = interp1(pdfs(ista).mm, pdfs(ista).pm, vsta(ista),...
        'linear', 0 ); % Probability of a velocity at a specific station, from our velocity model surface. 
    % Use linear interpolation for speed... Don't care for more precision
    % here I think .
    % Use 0 as extrap value, since we might sample outside a pdf. That's
    % fine: we just assume 0 probability. 
end

penalty = - sum(pv_mod); % Lower penalty is higher probability. 
penalty = penalty + roughness; 

end