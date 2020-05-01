function [fig, qs] = PlotT(T, fig)

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

% Retrieve the rotational and translational components of the transformation matrix.
[R, P] = TransToRp(T);

% Plot the given transformation matrix.
q1 = quiver3(P(1), P(2), P(3), R(1, 1), R(2, 1), R(3, 1), 'Linewidth', 2);
q2 = quiver3(P(1), P(2), P(3), R(1, 2), R(2, 2), R(3, 2), 'Color', q1.Color, 'Linewidth', 2);
q3 = quiver3(P(1), P(2), P(3), R(1, 3), R(2, 3), R(3, 3), 'Color', q1.Color, 'Linewidth', 2);

% Store the transformation matrix figure elements into an array.
qs = [q1 q2 q3];

end

