function Sigma = GetStressMatrix2D(sigma_xx, tau_xy, tau_yx, sigma_yy)

% This function returns the stress matrix associated with the given stress components.

% Inputs:
    % sigma_xx, tau_xy, tau_yx, sigma_yy = Scalar stress components of the stress matrix.

% Outputs:
    % Sigma = 2 x 2 matrix of the stress state of the form [sigma_xx tau_xy; tau_yx sigma_yy];.

% Retrieve the stress components from the stress matrix.
Sigma(1, 1) = sigma_xx;
Sigma(1, 2) = tau_xy;

Sigma(2, 1) = tau_yx;
Sigma(2, 2) = sigma_yy;

end

