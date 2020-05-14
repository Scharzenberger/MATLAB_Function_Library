function sigma_xx = GetTwoPlaneNormalBendingStress(My, Mz, y, z, Iy, Iz)

% This function computes the normal bending stress in a beam with applied moments My, Mz, second moment of the areas Iy, Iz, and cross-section position y, z.

% Note: x is along the length of the beam, y is along the height of the beam, z is along the width of the beam.

% Inputs:
    % My, Mz = Applied Moments about the y- and z-axis.
    % y, z = Cross-Section Position.
    % Iy, Iz = Second Moment of the Area about the y- and z-axes.

% Outputs:
    % sigma_xx = Normal Bending Stress on the x-plane in the x-direction.
    
% Compute the normal bending stress on the x-plane in the x-direction.
sigma_xx = -(Mz.*y)./Iz + (My.*z)./Iy;
    
end

