function [x_intersect] = funct_find_intersect(x, y, inters, options);
    arguments
        x
        y
        inters
        options.interpolate = true; 
    end 
% Find x values where y intersects inters.
% Can find multiple instances of intersection. 
% Does not need to be monotonic in y. 
% Can upsample using cubic interpolation. 
% Quick and dirty. 
% brb2023.04.03

if options.interpolate; 
    xold = x; 
    yold = y; 
    x = linspace(min(x), max(x), 10000); 
    y = interp1(xold, yold, x, "cubic"); 
end

inters_dist = y-inters; 
pos = inters_dist >= 0; 
sign_flip = pos(2:end) ~= pos(1:end-1); 

sign_flip = find(sign_flip); 

x_intersect = (x(sign_flip+1) + x(sign_flip))./2; 

end