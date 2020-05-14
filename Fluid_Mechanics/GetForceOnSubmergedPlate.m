function [F, x_CP, y_CP] = GetForceOnSubmergedPlate(p_CG, A, theta, Ixx, Ixy, gamma)

% This function computes the force F a static fluid applies on a plate of area A submerged at an angle theta, as well as the location at which this force acts.

% Inputs:
    % p_CG = Fluid Pressure at the Center of Gravity of the Plate
    % A = Area of the Plate
    % theta = Angle of the Plate with Respect to the Surface.
    % Ixx = Moment of Inertia of the Plate on the x-plane in the x-direction.
    % Ixy = Moment of Inertia of the Plate on the x-plane in the y-direction.
    % gamma = Specific Weight of the Fluid.
    
% Outputs:
    % F = Applied Force.
    % x_CP, y_CP = Center of Pressure Location with Respect to the Center of Gravity of the Plate.

% Compute the applied force.
F = p_CG.*A;

% Compute the center of pressure location.
x_CP = -gamma.*sin(theta).*(Ixy./F);
y_CP = -gamma.*sin(theta).*(Ixx./F);

end

