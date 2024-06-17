function z_smth = gaussian_smoothing_lab(ix,iy,longrid,latgrid,zgrid,wtgrid,smthdist,nstdinclude)
% z_smth = gaussian_smoothing_lab(ix,iy,longrid,latgrid,zgrid,wtgrid,smthdist,nstdinclude)
% function to apply a gaussian smoothing function across a surface,
% smoothing point ix, iy (i.e. grid(ix,iy)) using the points around it.
% smthdist is the sigma for the gaussian (i.e. 1/4 of the 95% interval.
% nstdinclude is the numnber of standard deviations *beyond which* the
% values are not incldued in the weighted sum (would have tiny weights
% anyway)

if nargin < 8
    nstdinclude = 2;
end

if isscalar(wtgrid)
    wtgrid = wtgrid*ones(size(longrid));
end

d2k = 111.1949;
% find distance IN KM to all other points and use to weight
dX = abs(longrid - longrid(ix,iy))*d2k*cosd(latgrid(ix,iy)); % account for sphericity in deg2km
dY = abs(latgrid - latgrid(ix,iy))*d2k;
dR2 = dX.^2 + dY.^2; % square distance, in km^2
dR2(dR2 > (nstdinclude.^2 * smthdist.^2) ) = nan; % kill distances larger than 2std
dRwt = exp(-dR2./(2*smthdist.^2));

nnan = ~isnan(zgrid) & ~isnan(wtgrid) & ~isnan(dRwt);

% calculate the smoothed value
z_smth = sum(zgrid(nnan).*dRwt(nnan).*wtgrid(nnan))./sum(dRwt(nnan).*wtgrid(nnan));
end


