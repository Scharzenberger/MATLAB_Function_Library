function tau_xy = GetTorsionalShearStress(T, r, J)

% This function computes the shear stress in a bar due to an applied torsional moment.

% Inputs:
    % T = Applied Torsional Moment.
    % r = Radial Position of Interest.
    % J = Polar Second Moment of the Area.
    
% Outputs:
    % tau_xy, tau_xz = Shear Stress on the x-plane in the y-direction (or z-direction, depending on the theta associated with the point of interest).
    
% Compute the angle of twist.
tau_xy = (T.*r)./J;

end

