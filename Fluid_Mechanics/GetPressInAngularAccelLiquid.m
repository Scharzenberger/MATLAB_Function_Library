function p = GetPressInAngularAccelLiquid(z, r, p0, g, rho, omega)

% This function computes the pressure p at depth z and radius r in a uniformly accelerating rotating liquid with density rho and constant angular velocity omega.

% Inputs:
    % p0 = Pressure in the Center of the Free Surface of the Rotating Fluid
    % z = Depth of Interest.
    % r = Radius of Interest.
    % g = Gravitational Constant.
    % rho = Density
    % omega = Constant Angular Velocity.

% Outputs:
    % p = Pressure at Point (z, r).

% Compute the pressure at the desired point.
p = p0 - rho.*g.*z + (1/2)*rho.*(r.^2).*(omega.^2);

end

