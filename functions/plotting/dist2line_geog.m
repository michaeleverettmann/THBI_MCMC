function [distance_perp,distance_par] = dist2line_geog( Q1,Q2,P,dr )
% [distance_perp,distance_par] = dist2line_geog( Q1,Q2,P )
%   calculates the distance perpendicular to and parallel to (as measured
%   from Q1) a line between geographic points Q1 [lat,lon] and Q2 [lat,lon]
%   of a series of points, P is a Nx2 vector of [lats,lons]. The
%   discretization along the line is dr (default is 0.01 deg).

if nargin < 4 || isempty(dr)
    dr = 0.01;
end

Np = size(P,1);

% % brute force: make lots of points between Q1 and Q2
% [qgc,qaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));
% rr = unique([0:dr:qgc,qgc])';
% Nq = length(rr);
% 
% QQ = zeros(Nq,2);
% for iq = 1:Nq
%     [QQ(iq,1), QQ(iq,2)] = reckon(Q1(1),Q1(2),rr(iq), qaz);
% end
% 
% DD = zeros(Np,Nq);
% for ip = 1:Np
% for iq = 1:Nq
%     DD(ip,iq) = distance(QQ(iq,1),QQ(iq,2),P(ip,1),P(ip,2));
% end
% end
% 
% [distance_perp,imin] = min(DD,[],2);
% if any(imin == 1) || any(imin == Nq)
%     warning('some points not between ends, but outside bounds')
% end
% distance_par = rr(imin);

%% cleverer method - use spherical trig...

% compute angle between line ends
[qgc,azq] = distance(Q1(1),Q1(2),Q2(1),Q2(2));

% Now cycle through points
%  work out az and distance from Q1 to each point
%  use sph trig to calc distance to line intersection (both perp and par)

distance_par = zeros(Np,1);
distance_perp = zeros(Np,1);

for ip = 1:Np
    [gci,azi] = distance(Q1(1),Q1(2),P(ip,1),P(ip,2));
    % relative azimuth to Pi from Q1-Q2 line
    theta = r2d(mp2pp(d2r(azq - azi)));
    % perp distance
    x1 = asind(sind(gci)*sind(theta));
    % par distance
    temp = sqrt( (sind(gci).^2 - 1)/(sind(gci).^2 * sind(theta).^2 -1) );
    x2 = min(abs([acosd(temp),acosd(-temp)]));
    % assign pos or neg
    if theta < 180 && theta > 0
        distance_perp(ip) = x1;
    else
        distance_perp(ip) = -x1;
    end
    if abs(theta) < 90
        distance_par(ip) = x2;
    else
        distance_par(ip) = -x2;
    end
end
    
    
    

    
    
    

end