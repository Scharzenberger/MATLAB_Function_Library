function tau_xy = GetShearBendingStress(V, Q, I, b)

% This function computes the shear bending stress in a beam with applied shear force V, first moment of the area Q, second moment of the area I, and beam width b.

% Note: x is along the length of the beam, y is along the height of the beam, z is along the width of the beam.

% Inputs:
    % M = Applied Moment.
    % Q = First Moment of the Area
    % I = Second Moment of the Area about the z-axis.
    % b = Width of Beam.

% Outputs:
    % tau_xy = Shear Bending Stress on the x-plane in the y-direction.
    
% Compute the shear bending stress on the x-plane in the y-direction.
tau_xy = (V.*Q)./(I.*b);
    
end

