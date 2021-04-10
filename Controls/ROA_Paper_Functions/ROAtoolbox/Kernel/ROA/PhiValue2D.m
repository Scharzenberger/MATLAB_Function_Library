function [meshxs, meshPhi] = PhiValue2D(projectoptions, fdata)
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
[meshxs, meshPhi, U, V] = gridnd2mesh(g, fdata, VF{:});

% Create a figure to store the plot.
figure('Color', 'w'), hold on, colormap('summer'), colorbar, xlabel('x_1'); ylabel('x_2'); zlabel('\phi(x, t)'); title('\phi(x, t) vs x_1, x_2')
% figure('Color', 'w'), hold on, colormap('jet'), colorbar, xlabel('x_1'); ylabel('x_2'); zlabel('\phi(x, t)'); title('\phi(x, t) vs x_1, x_2')

% Plot the phi surface.
surf(meshxs{:}, meshPhi, 'EdgeColor', 'none');

% Plot a 3D contour plot.
contour3(meshxs{:}, meshPhi, [0, 0], 'r', 'Linewidth', 2);

%Plot the associated flow lines.
h = streamslice(meshxs{:}, U, V, 0.4);

% Set the flow plot color.
set(h, 'color', 'b')
% set(h, 'color', 'k')

% Plot the flow outside of the ROA.
for i=1:length(h)
    Xq = get(h(i), 'xdata');
    Yq = get(h(i), 'ydata');
    Vq = interp2(meshxs{:}, meshPhi, Xq, Yq);
    set(h(i),'zdata', Vq);
end


