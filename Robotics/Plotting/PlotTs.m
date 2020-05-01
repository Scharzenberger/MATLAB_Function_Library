function [fig, qs] = PlotTs(Ts, fig)

% Set the default arguments.
if nargin < 2, fig = []; end

% Determine whether to create a new figure or select an existing one.
if isempty(fig)             % If no figure was provided...
    
    % Create a new figure.
    fig = figure('color', 'w'); hold on
    
else                        % Otherwise...
    
    % Select the provided figure.
    figure(fig)
    
end

% Get the color order associated with this figure.
colorOrder = get(gca, 'ColorOrder');

% Retrieve the number of childern on the plot.
num_childern = length( get(gca, 'Children') );

% Retrieve the previous color index.
color_index = mod( num_childern - 1, size(colorOrder, 1) ) + 1;

% Define the dimension of this problem.
dim = 3;                    % [#] Number of dimensions.  This is fixed due to the nature of the transformation matrices.

% Retrieve size information from the transformation matrices.
num_joints = size(Ts, 3);                   % [#] Number of Joints.
num_angles = size(Ts, 4);                   % [#] Number of Angles.

% Create a variable to store the transformation matrix figure elements.
qs = gobjects(dim, num_joints, num_angles);

% Creat an array to store the used colors.
used_colors = zeros(num_joints, 3);

% Create a counter variable.
counter = color_index;

% Plot each of the transformation matrices.
for k1 = 1:num_angles                       % Iterate through each of the angles...
    for k2 = 1:num_joints                   % Iterate through each of the joints...
        
        % Determine whether to select a new color or use an existing color.
        if k1 == 1                      % If this is the first angle...
            
            % Advance the counter.
            counter = counter + 1;
            
            % Compute the color index.
            color_index = mod( counter - 1, size(colorOrder, 1) ) + 1;
            
            % Retrieve the current color.
            color = colorOrder(color_index, :);
            
            % Store the used colors.
            used_colors(k2, :) = color;
            
        else                            % Otherwise...
            
            % Use one of the previously used colors.
            color = used_colors(k2, :);
            
        end
        
        % Retrieve the rotational and translational components of the transformation matrix.
        [R, P] = TransToRp(Ts(:, :, k2, k1));
        
        % Create an arrow for each quiver.
        for k3 = 1:dim                  % Iterate through each dimension...
            
            % Plot the given transformation matrix.
            qs(k3, k2, k1) = quiver3(P(1), P(2), P(3), R(1, k3), R(2, k3), R(3, k3), 'Color', color, 'Linewidth', 2);
            
        end
        
    end
end

end

