function Keq = GetCantileveredBeamEqStiffness(E, Iarea, L)

% This function computes the equivalent stiffness of a cantilevered beam for displacement at the end of the beam.

% Assumptions:
    % This is only valid for small beam tip deflections.

% Inputs:
    % E = Young's Modulus of Elasticity
    % Iarea = Second Moment of the Area.
    % L = Length of the Beam.
    
% Outputs:
    % Keq = Equivalent Spring Stiffness.

% Compute the equivalent spring stiffness.
Keq = 3*E.*Iarea./(L.^3);
    
end

