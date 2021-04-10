function [h, VecField] = PhasePortrait(ProjectOptions, density)
% Compute the phase portrait and draw it.
% Since toolboxLS-1.1.1 uses ndgrid, one should transfer the grid to
%   mexhgrid-based array.
% Parameters:
%   projectoptions  Project Options for ROAtoolbox
%   densit    the density of trajectories
%   VecField        cell contains the stream.
%   h    handle of trajectories
%
% YUAN Guoqiang, Oct, 2016
%

% Retrieve the grid from the project options.
g = ProjectOptions.Grid;

% Retrieve the system vector field from the project options.
VF = ProjectOptions.VectorField;

% Retrieve the ROA initial condition.
sep = ProjectOptions.InitCircleCen;

% Determine how to handle default input arguments.
if nargin < 2               % If there are less than two inputs...
    switch g.dim
        case 2
            density = 2;
        case 3
            density = 0.4;
        case 4
            density = 2;
        case 5
            density = 2;
    end
end

% Determine how to visualize the phase portrait.
switch (g.dim)
    case 2
        
        % Convert the nd grid to a mesh grid.
        [meshxs, U, V] = gridnd2mesh(g, VF{1}, VF{2});
        
        % Retrieve the vector field.
        VecField = {U; V};
        
        % Retrieve the grid meshes.
        MX = meshxs{1}; MY = meshxs{2};
        
        % Create a figure to store the phase portrait.
        figure('Color', 'w'), hold on, box on, xlabel('x_1'), ylabel('x_2'), title('\phi(x, t) vs x_1, x_2'), axis equal, axis(g.axis)
        
        % Plot the 2D phase portrait.
        h = streamslice(MX, MY, U, V, density);
        
    case 3
        
        % Convert the nd grid to a mesh grid.
        [meshxs, U, V, W] = gridnd2mesh(g, VF{1}, VF{2}, VF{3});
        
        % Store the meshgrid vector field.
        VecField = {U; V; W};
        
        % Retrieve the grid meshes.
        MX = meshxs{1}; MY = meshxs{2}; MZ = meshxs{3};
        figure('Color', 'w'), hold on, box on, xlabel('x_1'), ylabel('x_2'), zlabel('x_3'), axis equal, axis(g.axis)
        
        % Plot the streamslices of the phase portrait.
        h = streamslice(MX, MY, MZ, U, V, W, sep(1), [], [], density);
        streamslice(MX, MY, MZ, U, V, W, [], sep(2), [], density);
        streamslice(MX, MY, MZ, U, V, W, [], [], sep(3), density);
        
    case 4
        
        % Convert the nd grid to a mesh grid.
        [meshxs, U1, V1, U2, V2] = gridnd2mesh(g, VF{1}, VF{2}, VF{3}, VF{4});
        
        % Retrieve the vector field.
        VecField = {U1; V1; U2; V2};
        
        % Retrieve the grid meshes.
        MX1 = meshxs{1}; MY1 = meshxs{2}; MX2 = meshxs{3}; MY2 = meshxs{4};

        % Project the grid meshes to appropriate lower dimensions.
        MX1 = MX1(:, :, round(size(MX1, 3)/2), round(size(MX1, 4)/2));
        MY1 = MY1(:, :, round(size(MY1, 3)/2), round(size(MY1, 4)/2));
        U1 = U1(:, :, round(size(U1, 3)/2), round(size(U1, 4)/2));
        V1 = V1(:, :, round(size(V1, 3)/2), round(size(V1, 4)/2));

        MX2 = reshape( MX2(round(size(MX2, 1)/2), round(size(MX2, 2)/2), :, :), [size(MX2, 3), size(MX2, 4)] )';
        MY2 = reshape( MY2(round(size(MY2, 1)/2), round(size(MY2, 2)/2), :, :), [size(MY2, 3), size(MY2, 4)] )';
        U2 = reshape( U2(round(size(U2, 1)/2), round(size(U2, 2)/2), :, :), [size(U2, 3), size(U2, 4)] )';
        V2 = reshape( V2(round(size(V2, 1)/2), round(size(V2, 2)/2), :, :), [size(V2, 3), size(V2, 4)] )';
        
        % Create a figure to store the phase portrait.
        h1 = figure('Color', 'w'); hold on, box on, xlabel('x_1'), ylabel('x_2'), title('Flow vs x1 & x2 (x3 = 0, x4 = 0)'), axis equal, axis([min(min(MX1)), max(max(MX1)), min(min(MY1)), max(max(MY1))]), streamslice(MX1, MY1, U1, V1, density);
        h2 = figure('Color', 'w'); hold on, box on, xlabel('x_3'), ylabel('x_4'), title('Flow vs x3 & x4 (x1 = 0, x2 = 0)'), axis equal, axis([min(min(MX2)), max(max(MX2)), min(min(MY2)), max(max(MY2))]), streamslice(MX2, MY2, U2, V2, density);

        % Store the plot handles.
        h = [h1, h2];
        
        case 5
        
        % Convert the nd grid to a mesh grid.
        [meshxs, U1, V1, U2, V2] = gridnd2mesh(g, VF{1}, VF{2}, VF{3}, VF{4});
        
        % Retrieve the vector field.
        VecField = {U1; V1; U2; V2};
        
        % Retrieve the grid meshes.
        MX1 = meshxs{1}; MY1 = meshxs{2}; MX2 = meshxs{3}; MY2 = meshxs{4};

        % Remove the extra dimension.
        MX1 = MX1(:, :, :, :, 1); MY1 = MY1(:, :, :, :, 1); MX2 = MX2(:, :, :, :, 1); MY2 = MY2(:, :, :, :, 1);
        U1 = U1(:, :, :, :, 1); V1 = V1(:, :, :, :, 1); U2 = U2(:, :, :, :, 1); V2 = V2(:, :, :, :, 1);
        
        % Project the grid meshes to appropriate lower dimensions.
        MX1 = MX1(:, :, round(size(MX1, 3)/2), round(size(MX1, 4)/2));
        MY1 = MY1(:, :, round(size(MY1, 3)/2), round(size(MY1, 4)/2));
        U1 = U1(:, :, round(size(U1, 3)/2), round(size(U1, 4)/2));
        V1 = V1(:, :, round(size(V1, 3)/2), round(size(V1, 4)/2));

        MX2 = reshape( MX2(round(size(MX2, 1)/2), round(size(MX2, 2)/2), :, :), [size(MX2, 3), size(MX2, 4)] )';
        MY2 = reshape( MY2(round(size(MY2, 1)/2), round(size(MY2, 2)/2), :, :), [size(MY2, 3), size(MY2, 4)] )';
        U2 = reshape( U2(round(size(U2, 1)/2), round(size(U2, 2)/2), :, :), [size(U2, 3), size(U2, 4)] )';
        V2 = reshape( V2(round(size(V2, 1)/2), round(size(V2, 2)/2), :, :), [size(V2, 3), size(V2, 4)] )';
        
        % Create a figure to store the phase portrait.
%         h1 = figure('Color', 'w'); hold on, box on, xlabel('x_1'), ylabel('x_2'), title('Flow vs x1 & x2 (x3 = 0, x4 = 0)'), axis equal, axis([min(min(MX1)), max(max(MX1)), min(min(MY1)), max(max(MY1))]), streamslice(MX1, MY1, U1, V1, density);
%         h2 = figure('Color', 'w'); hold on, box on, xlabel('x_3'), ylabel('x_4'), title('Flow vs x3 & x4 (x1 = 0, x2 = 0)'), axis equal, axis([min(min(MX2)), max(max(MX2)), min(min(MY2)), max(max(MY2))]), streamslice(MX2, MY2, U2, V2, density);

        h1 = figure('Color', 'w'); hold on, box on, xlabel('x_1'), ylabel('x_2'), title('Flow vs x1 & x2 (x3 = 0, x4 = 0)'), axis([min(min(MX1)), max(max(MX1)), min(min(MY1)), max(max(MY1))]), streamslice(MX1, MY1, U1, V1, density);
        h2 = figure('Color', 'w'); hold on, box on, xlabel('x_3'), ylabel('x_4'), title('Flow vs x3 & x4 (x1 = 0, x2 = 0)'), axis([min(min(MX2)), max(max(MX2)), min(min(MY2)), max(max(MY2))]), streamslice(MX2, MY2, U2, V2, density);


        % Store the plot handles.
        h = [h1, h2];
        
    otherwise
        
        % Throw an error.
        error('Can not draw phase portrait for system with dimention: %s!', g.dim);
        
end