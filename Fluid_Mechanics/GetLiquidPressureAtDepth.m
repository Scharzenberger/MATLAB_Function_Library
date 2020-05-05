function p2 = GetLiquidPressureAtDepth(p1, z1, z2, gamma)

% This function computes the pressure at depth z2 given the pressure p1 at depth z1 and the specific weight of the liquid of interest.

% Assumptions:
    % The density of the liquid does not vary significantly with depth.  Since liquids are nearly incompressible, this is more or less valid for liquids.

% Inputs:
    % p1 = Pressure at Depth z1.
    % z1, z2 = Depths of Points 1 and 2.
    % gamma = Specific Weight of the Fluid.

% Outputs:
    % p2 = Pressure at Depth z2.

% Compute the pressure at the given depth.
p2 = p1 - gamma*(z2 - z1);
    
end

