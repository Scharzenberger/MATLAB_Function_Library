function gamma = Density2SpecificWeight(rho, g)

% This function computes the specific weight associated with a given density rho.

% Inputs:
    % rho = Fluid Density.
    % g = Gravitational Constant.

% Outputs:
    % gamma = Specific Weight.
    
% Set the default input arguments.
if nargin < 2, g = 9.81; end
 
% Compute the specific weight.
gamma = rho.*g;

end

