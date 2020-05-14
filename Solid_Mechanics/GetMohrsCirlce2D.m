function [C, R] = GetMohrsCirlce2D(Sigma, bMakePlot)

% This function computes the center location and radius of the 2D Mohr's Circle (i.e., Mohr's Circle for Plane Stress) associated with a stress state Sigma.

% Inputs:
    % Sigma = 2 x 2 matrix of the form [sigma_x; tau_xy; tau_yx sigma_y]; that describes the plane stress state at the point of interest.
    % bMakePlot = Boolean that determines whether the Mohr's Circle for Plane Stress is plotted.

% Outputs:
    % C = 2 x 1 array of the form [sigma; tau] that describes the center of Mohr's Circle for Plane Stress.
    % R = Scalar of the form sigma that describes the radius of Mohr's Circle for Plane Stress.
   
% Ensure that the stress state is of the appropriate size and is symmetric.
if any(size(Sigma) ~= [2 2]), error('Sigma must be a 2 x 2 matrix describing a state of plane stress.\n'); end
if Sigma(1, 2) ~= Sigma(2, 1), error('Sigma must be symmetric.\n'); end

% Set the default input arguments.
if nargin < 2, bMakePlot = true; end

% Retrieve the components of the stress matrix.
[sigma_x, tau_xy, ~, sigma_y] = GetStressComponents2D(Sigma);
    
% Compute the center of Mohr's Circle for Plane Stress.
C = [(sigma_x + sigma_y)/2; 0];

% Compute the radius of Mohr's Circle for Plane Stress.
R = sqrt( ( (sigma_x - sigma_y)/2 ).^2 + tau_xy.^2 );
    
% Determine whether to plot Mohr's Circle for Plane Stress.
if bMakePlot                % If we are asked to plot Mohr's Circle for Plane Stress...

    % Create a parameter for Mohr's Cirlce for Plane Stress.
    ts = linspace(0, 2*pi, 100);
    
    % Compute the points on Mohr's Circle for Plane Stress.
    xs = R*cos(ts) + C(1);
    ys = R*sin(ts) + C(2);
    
    % Create a figure for Mohr's Circle for Plane Stress.
    figure('Color', 'w', 'Name', 'Mohr''s Circle for Plane Stress'), hold on, grid on, xlabel('Normal Stress, $\sigma$ [Pa]', 'Interpreter', 'Latex'), ylabel('Shear Stress, $\tau$ [Pa]', 'Interpreter', 'Latex'), title('Mohr''s Circle for Plane Stress')
    
    % Plot Mohr's Circle for Plane Stress.
    plot(xs, ys, '-', 'Linewidth', 3)
    
end

end

