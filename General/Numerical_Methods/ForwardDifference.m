function gradf = ForwardDifference(fs, grid, n)

% Retrieve the size of the specified grid.
grid_size = size(grid);

% Retrieve the number of dimensions.
num_dims = grid_size(end);

% Remove the last entry from the grid size variable (because this is the number of dimensions, not the grid size per dimension).
grid_size(end) = [];

% Define the total number of elements in the grid.
num_entries = numel(fs);

gradf = zeros([grid_size - n, num_dims]);

for k = 1:num_entries

    % Compute the gradient at this location.
    gradfk = FowardDifferenceAtPoint(fs, grid, n, k);

%     % Store the gradient at this location into the gradient array.
%     gradf(:, :, 
    
end

end

