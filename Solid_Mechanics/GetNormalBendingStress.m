function sigma_xx = GetNormalBendingStress(M, y, I)

% This function computes the normal bending stress in a beam with applied moment M, second moment of the area I, and cross-section height y.

% Note: x is along the length of the beam, y is along the height of the beam, z is along the width of the beam.

% Inputs:
    % M = Applied Moment.
    % y = Cross-Section Height.
    % I = Second Moment of the Area about the z-axis.

% Outputs:
    % sigma_xx = Normal Bending Stress on the x-plane in the x-direction.
    
% Compute the normal bending stress on the x-plane in the x-direction.
sigma_xx = -(M.*y)./I;
    
end

