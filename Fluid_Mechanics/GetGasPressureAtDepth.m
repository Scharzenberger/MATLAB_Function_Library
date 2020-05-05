function p2 = GetGasPressureAtDepth(p1, z1, z2, g, R, T0)

% This function computes the pressure at depth z2 given the pressure p1 at depth z1, the associated gas constant R, and the temperature of the gas.

% Assumptions:
    % The temperature of the gas is the same at the two depths.

% Inputs:
    % p1 = Pressure at Depth z1.
    % z1, z2 = Depths of Points 1 and 2.
    % g = Gravitational Constant
    % R = Gas Constant
    % T0 = Gas Temperature
    
% Outputs:
    % p2 = Pressure at Depth z2.

% Compute the pressure at the given depth.
p2 = p1*exp(-(g*(z2 - z1))/(R*T0));
    
end


