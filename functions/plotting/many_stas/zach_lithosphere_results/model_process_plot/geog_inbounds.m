function [inbounds] = geog_inbounds(lat,lon)
%[inbounds] = geog_inbounds(lat,lon)
%   Function to evaluate whether a point (lat,lon can be NxM matrices, or
%   scalars) is within the lat/lon bounds for which we have long-period
%   surface wave data. I.e. specifically not including most of Florida or
%   the narrow New England area.

szla = size(lat);
szlo = size(lon);

if ~isequal(szlo,szla)
    if numel(lat)~=1 && numel(lon)~=1 
        error('lat and lon must be same size, or >=one must be scalar')
    end
end

% make same size, if not
if numel(lat)==1 && numel(lon)~=1
    lat = lat* ones(size(lon));
elseif numel(lon)==1 && numel(lat)~=1
    lon = lon* ones(size(lat));
end

% set up. All in by default
inbounds = true(size(lat));

%kill FL
inbounds(lat<30.5) = false;
inbounds(lat<31 & lon>-82.2) = false;
%kill NE
inbounds(lat<42.25 & lon>-72) = false;
%kill bad ME
inbounds(lon>-68.5) = false;
% kill offshore ME/NH
[~,azc] = distance_deg_nomap(42.,-71.3,44.3,-68.7);
[~,az] = distance_deg_nomap(42.,-71.3,lat,lon);
inbounds( (lat > 42.) & (lat < 45.) & (az > azc) & (az < 120) ) = false;
% if  (az > azc) & (az < 120) 
%     keyboard
% end

end

