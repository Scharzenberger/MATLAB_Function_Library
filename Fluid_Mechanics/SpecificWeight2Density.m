function rho = SpecificWeight2Density(gamma, g)

% This function computes the density rho associated with a given specific weight gamma.

% Inputs:
    % gamma = Specific Weight.
    % g = Gravitational Constant.

% Outputs:
    % rho = Fluid Density.

% Set the default input arguments.
if nargin < 2, g = 9.81; end
 
% Compute the specific weight.
rho = gamma./g;

end

