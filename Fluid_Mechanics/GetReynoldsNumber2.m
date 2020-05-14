function Re = GetReynoldsNumber2(V, L, nu)

% This function computes the Reynold's Number associated with a flow given basic flow properties, among which is the fluid kinematic viscosity.

% Inputs:
    % V = Fluid Velocity
    % L = Characteristic Length (e.g., Pipe Diameter)
    % nu = Fluid Kinematic Viscosity

% Outputs:
    % Re = Reynold's Number
    
% Compute the Reynold's Number associated with the flow.
Re = V*L/nu;

end

