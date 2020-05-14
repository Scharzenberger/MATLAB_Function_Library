function F = GetBuoyantForce(rho, g, V)

% This function computes the buoyant force acting on a (partially) submerged body.

% Inputs:
    % rho = Fluid Density
    % g = Gravitational Constant
    % V = Volume of Displaced Fluid

% Outputs:
    % F = Buoyant Force
    
% Compute the buoyant force.
F = rho*g*V;

end

