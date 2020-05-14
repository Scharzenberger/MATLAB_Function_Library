function [p2, theta] = GetPressInLinearAccelLiquid(p1, z1, z2, rho, g, ax, az)

% This function computes the pressure p2 of a liquid of density rho at depth z2 when the fluid is accelerating uniformly at ax, az.

% Inputs:
    % p1 = Pressure at Depth z1.
    % z1, z2 = Depths 1 and 2.
    % rho = Liquid Density
    % g = Gravitational Constant
    % ax, az = Acceleration in the x & z directions.

% Outputs:
    % p2 = Pressure at Depth z2.
    % theta = Angle of the Free Surface.

% Compute the total uniform acceleration magnitude.
G = sqrt( ax.^2 + (g + az.^2) );
    
% Compute the angle of the free surface.
theta = atan2(ax, g + az);

% Compute the pressure at depth z2.
p2 = p1 + rho.*G.*(z2 - z1).*csc(theta);


end

