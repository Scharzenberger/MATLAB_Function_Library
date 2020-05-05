function Re = GetReynoldsNumber1(rho, V, L, mu)

% This function computes the Reynold's Number associated with a flow given basic flow properties, among which is the fluid viscosity.

% Inputs:
    % rho = Fluid Density
    % V = Fluid Velocity
    % L = Characteristic Length (e.g., Pipe Diameter)
    % mu = Fluid Viscosity

% Outputs:
    % Re = Reynold's Number
    
% Compute the Reynold's Number associated with the flow.
Re = rho*V*L/mu;

end

