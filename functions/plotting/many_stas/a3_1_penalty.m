function [penalty] = a3_1_penalty(vgrid)
global pdf rough_scale dx2 dy2 xgrid ygrid stax stay

%% 

% Use an example pdf while developing surface inversion stuff. 
% fpdf = '~/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/pdf_example.mat'; 
% pdf = load(fpdf).pdf_example; 
% rough_scale = 5e-11; 

vgrid(isnan(vgrid)) = nanmean(nanmean(vgrid)); 

% dx = xgrid(2:end,:)-xgrid(1:end-1,:); 
% dy = ygrid(:,2:end)-ygrid(:,1:end-1); 
% dx2 = xgrid(1:end-2,:) - 2*xgrid(2:end-1,:) + xgrid(3:end,:);
% dy2 = xgrid(:,1:end-2) - 2*ygrid(1,2:end-1) + ygrid(:,3:end); 
% dx2 = ((xgrid(1:end-2,:) - xgrid(3:end,:))/2).^2; 
% dy2 = ((ygrid(:,1:end-2) - ygrid(:,3:end))/2).^2; 
% x_dx2 = xgrid(2:end-1,:); y_dx2 = ygrid(2:end-1,:); 
% y_dy2 = ygrid(:,2:end-1); x_dy2 = xgrid(:,2:end-1); 

% Options for optimizing. 
% Use one value for dx and dy. But as is, it's more versatile for map
% projections. 


dvdx2 = (vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2; 
dvdy2 = (vgrid(:,1:end-2) - 2*vgrid(:,2:end-1) + vgrid(:,3:end))./dy2; 
% dvdx2(isnan(dvdx2)) = 2; 
% dvdy2(isnan(dvdy2)) = 2; 


% figure(2); clf; hold on; 
% contourf(x_dx2, y_dx2, dvdx2, 50); 
% colorbar(); 

roughness = sum(dvdx2 .^ 2, 'all') + sum(dvdy2 .^ 2, 'all'); 
roughness = roughness * rough_scale; 

vsta = interpn(xgrid, ygrid, vgrid, stax, stay, 'linear'); % Modelled velocity at each station 

% ista = 1; 
pv_mod = zeros(size(stax)); 
for ista=1:length(pv_mod); 
    pv_mod(ista) = interp1(pdf.v, pdf.p, vsta(ista) ); % Probability of a velocity at a specific station, from our velocity model surface. 
end

penalty = - sum(pv_mod); % Lower penalty is higher probability. 
penalty = penalty + roughness; 

% pdf.p

% fun(vgrid)@(); 

% global 

end