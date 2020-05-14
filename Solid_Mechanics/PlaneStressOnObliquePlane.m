function [sigma_normal, tau_tangent] = PlaneStressOnObliquePlane(Sigma, phi)

% This function computes the normal and shear stresses acting on an oblique plane whose normal is at an angle of phi with respect to the x direction at a point with plane stress state Sigma.

% Inputs:
    % Sigma = 2 x 2 matrix of the form [sigma_x; tau_xy; tau_yx sigma_y]; that describes the plane stress state at the point of interest.
    % phi = Angle of plane on which to compute the stress state.
    
% Outputs:
    % sigma_normal = Normal stress component normal to the plane.
    % tau_tangent = Shear stress component tangent to the plane.
    
% Ensure that the stress state is of the appropriate size and is symmetric.
if any(size(Sigma) ~= [2 2]), error('Sigma must be a 2 x 2 matrix describing a state of plane stress.\n'); end
if Sigma(1, 2) ~= Sigma(2, 1), error('Sigma must be symmetric.\n'); end

% Retrieve the stress components.
[sigma_x, tau_xy, ~, sigma_y] = GetStressComponents2D(Sigma);
    
% Compute the normal stress.
sigma_normal = (sigma_x + sigma_y)/2 + ((sigma_x - sigma_y)/2)*cos(2*phi) + tau_xy*sin(2*phi);

% Compute the shear stress.
tau_tangent = -((sigma_x - sigma_y)/2)*sin(2*phi) + tau_xy*cos(2*phi);


end

