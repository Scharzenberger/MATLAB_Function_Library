function [sigma_ps, tau_ps] = GetMohrsCircle3D( Sigma, bMakePlot )

% This function computes the principal normal stress and principal shear stresses associated with a given stress state and plots the associated 3D Mohr's Cirlce

% Inputs:
    % Sigma = 3 x 3 matrix of the form [sigma_xx tau_xy tau_xz; tau_yx sigma_yy tau_yz; tau_zx tau_zy sigma_zz]; that describes the stress state at the point of interest.
    % bMakePlot = Boolean that determines whether the Mohr's Circle for Plane Stress is plotted.

% Outputs:
    % sigma_ps = 3 x 1 column vector of principal normal stresses.
    % tau_ps = 3 x 1 column vector of principal shear stresses.
    
% Ensure that the stress state is of the appropriate size and is symmetric.
if any(size(Sigma) ~= [3 3]), error('Sigma must be a 3 x 3 matrix describing a state of plane stress.\n'); end
if any(Sigma ~= Sigma'), error('Sigma must be symmetric.\n'); end

% Set the default input arguments.
if nargin < 2, bMakePlot = true; end

% Compute the principal stresses.
sigma_ps = eig(Sigma);
    
% Compute the principal shear stresses.
tau_12 = (sigma_ps(1) - sigma_ps(2))/2;
tau_23 = (sigma_ps(2) - sigma_ps(3))/2;
tau_13 = (sigma_ps(1) - sigma_ps(3))/2;

% Store the principal shear stresses as a column vector.
tau_ps = [tau_12; tau_23; tau_13];

% Determine whether to plot the 3D Mohr's Cirlce.
if bMakePlot                % If the user requests that we plot the 3D Mohr's Circle...
    
    % Create a parameter for the circles.
    ts = linspace(0, 2*pi, 100);
    
    % Define the radii of the three circles.
    R12 = tau_12;
    R23  = tau_23;
    R13 = tau_13;
    
    % Compute the centers of the three circles.
    C12 = sigma_ps(2) + R12;
    C23 = sigma_ps(3) + R23;
    C13 = sigma_ps(3) + R13;
    
    % Create the points for the three cirlces.
    xs12 = R12*cos(ts) + C12; ys12 = R12*sin(ts);
    xs23 = R23*cos(ts) + C23; ys23 = R23*sin(ts);
    xs13 = R13*cos(ts) + C13; ys13 = R13*sin(ts);

    % Create a figure for the 3D Mohr's Circle plot.
    figure('Color', 'w', 'Name', '3D Mohr''s Cirlce'), hold on, grid on, xlabel('Normal Stress, $\sigma$ [Pa]', 'Interpreter', 'Latex'), ylabel('Shear Stress, $\tau$ [Pa]', 'Interpreter', 'Latex'), title('3D Mohr''s Circle')
    
    % Plot the three circles.
    plot(xs12, ys12, '-', 'Linewidth', 3)
    plot(xs23, ys23, '-', 'Linewidth', 3)
    plot(xs13, ys13, '-', 'Linewidth', 3)

    % Plot the prinipcal stresses.
    plot(sigma_ps, zeros(size(sigma_ps)), '.', 'Markersize', 20)
    
end

end

