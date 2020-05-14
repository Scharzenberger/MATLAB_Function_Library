function Sigma = GetStressMatrix3D(sigma_xx, tau_xy, tau_xz, tau_yx, sigma_yy, tau_yz, tau_zx, tau_zy, sigma_zz)

% This function returns the stress matrix associated with the given stress components.

% Inputs:
    % sigma_xx, tau_xy, tau_xz, tau_yx, sigma_yy, tau_yz, tau_zx, tau_zy, sigma_zz = Scalar stress components of the stress matrix.

% Outputs:
    % Sigma = 3 x 3 matrix of the stress state of the form [sigma_xx tau_xy tau_xz; tau_yx sigma_yy tau_yz; tau_zx tau_zy sigma_zz];.
    
% Retrieve the stress components from the stress matrix.
Sigma(1, 1) = sigma_xx;
Sigma(1, 2) = tau_xy;
Sigma(1, 3) = tau_xz;

Sigma(2, 1) = tau_yx;
Sigma(2, 2) = sigma_yy;
Sigma(2, 3) = tau_yz;

Sigma(3, 1) = tau_zx;
Sigma(3, 2) = tau_zy;
Sigma(3, 3) = sigma_zz;

end

