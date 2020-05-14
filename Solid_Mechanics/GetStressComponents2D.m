function [sigma_xx, tau_xy, tau_yx, sigma_yy] = GetStressComponents2D(Sigma)

% This function retrieve the individual stress components of a stress matrix.

% Inputs:
    % Sigma = 2 x 2 matrix of the stress state of the form [sigma_xx tau_xy; tau_yx sigma_yy];.

% Outputs:
    % sigma_xx, tau_xy, tau_yx, sigma_yy = Scalar stress components of the stress matrix.
    
% Retrieve the stress components from the stress matrix.
sigma_xx = Sigma(1, 1);
tau_xy = Sigma(1, 2);

tau_yx = Sigma(2, 1);
sigma_yy = Sigma(2, 2);

end

