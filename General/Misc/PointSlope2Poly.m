function poly = PointSlope2Poly(point, slope)

% This function computes the linear polynomial that passes through the point "point" with slope "slope."

% Compute the linear polynomial that passes through the point "point" with slope "slope."
poly = [slope,  point(2) - slope*point(1)];


end

