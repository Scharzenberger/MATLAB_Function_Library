function [meshxs, meshPhi] = PhiValue4D(projectoptions, fdata)
% Draw the level set function phi.
% Since toolboxLS-1.1.1 uses ndgrid, one should transfer the grid to
%   mexhgrid-based array.
% Parameters:
%   projectoptions  Project Options for ROAtoolbox
%   fdata        array contains the function value on the grid.
%
% YUAN Guoqiang, Oct, 2016
%

% Retrieve the grid.
g = projectoptions.Grid;

% Retrieve the flow field.
VF = projectoptions.VectorField;

% Convert the nd grid to a mesh grid.
[meshxs, meshPhi, U1, V1, U2, V2] = gridnd2mesh(g, fdata, VF{:});

% Retrieve the girds.
X1 = meshxs{1}; X2 = meshxs{2}; X3 = meshxs{3}; X4 = meshxs{4};

X1_proj = X1(:, :, round(size(X1, 3)/2), round(size(X1, 4)/2));
X2_proj = X2(:, :, round(size(X2, 3)/2), round(size(X2, 4)/2));
U1_proj = U1(:, :, round(size(U1, 3)/2), round(size(U1, 4)/2));
V1_proj = V1(:, :, round(size(V1, 3)/2), round(size(V1, 4)/2));
Phi1_proj = meshPhi(:, :, round(size(meshPhi, 3)/2), round(size(meshPhi, 4)/2));

X3_proj = reshape( X3(round(size(X3, 1)/2), round(size(X3, 2)/2), :, :), [size(X3, 3), size(X3, 4)] )';
X4_proj = reshape( X4(round(size(X4, 1)/2), round(size(X4, 2)/2), :, :), [size(X4, 3), size(X4, 4)] )';
U2_proj = reshape( U2(round(size(U2, 1)/2), round(size(U2, 2)/2), :, :), [size(U2, 3), size(U2, 4)] )';
V2_proj = reshape( V2(round(size(V2, 1)/2), round(size(V2, 2)/2), :, :), [size(V2, 3), size(V2, 4)] )';
Phi2_proj = reshape( meshPhi(round(size(meshPhi, 1)/2), round(size(meshPhi, 2)/2), :, :), [size(meshPhi, 3), size(meshPhi, 4)] )';

% Create a figure to store the plot.
figure('Color', 'w'), hold on, colormap('summer'), colorbar, xlabel('x_3'); ylabel('x_4'); zlabel('\phi(x,t)'); title('\phi(x,t) vs x_3, x_4 (x_1 = 0, x_2 = 0)')

% Plot the phi surface.
surf(X3_proj, X4_proj, Phi2_proj, 'EdgeColor', 'none');

% Plot a 3D contour plot.

contour3(X3_proj, X4_proj, Phi2_proj, [1, 1], 'r', 'Linewidth', 2);

%Plot the associated flow lines.
h = streamslice(X3_proj, X4_proj, U2_proj, V2_proj, 0.4);

% Set the flow plot color.
set(h, 'color', 'b')

% Correct the z value of the flow lines.
for i=1:length(h)
    
    % Retrieve the x & y data associated with this flow line.
    xs_line = get(h(i), 'xdata'); ys_line = get(h(i), 'ydata');
    
    % Compute the z data for this flow line.
    zs_line = interp2(X3_proj, X4_proj, Phi2_proj, xs_line, ys_line);
    
    % Update the z data for this flow line.
    set(h(i),'zdata', zs_line);
    
end


