function [sigma_xx, tau_xy, tau_xz, tau_yx, sigma_yy, tau_yz, tau_zx, tau_zy, sigma_zz] = GetStressComponents3D(Sigma)

% This function retrieve the individual stress components of a stress matrix.

% Inputs:
    % Sigma = 3 x 3 matrix of the stress state of the form [sigma_xx tau_xy tau_xz; tau_yx sigma_yy tau_yz; tau_zx tau_zy sigma_zz];.

% Outputs:
    % sigma_xx, tau_xy, tau_xz, tau_yx, sigma_yy, tau_yz, tau_zx, tau_zy, sigma_zz = Scalar stress components of the stress matrix.
    
% Retrieve the stress components from the stress matrix.
sigma_xx = Sigma(1, 1);
tau_xy = Sigma(1, 2);
tau_xz = Sigma(1, 3);

tau_yx = Sigma(2, 1);
sigma_yy = Sigma(2, 2);
tau_yz = Sigma(2, 3);

tau_zx = Sigma(3, 1);
tau_zy = Sigma(3, 2);
sigma_zz = Sigma(3, 3);

end

