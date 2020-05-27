function gradf = ForwardDifference(fs, grid, n)

% Retrieve the size of the specified grid.
grid_size = size(grid);

% Retrieve the number of dimensions.
num_dims = grid_size(end);

% Remove the last entry from the grid size variable (because this is the number of dimensions, not the grid size per dimension).
grid_size(end) = [];

% Define the total number of elements in the grid.
num_entries = numel(fs);

% Initialize the size of the gradient array.
gradf = zeros([grid_size - n, num_dims]);

% Compute the gradient at each point in the grid.
for k = 1:num_entries               % Iterate through each point in the grid...
    
    % Retrieve the dimensional indexes associated with this entry.
    k_indexes = AbsInd2DimInd( grid_size, k );
        
    % Determine whether we can perform a forward difference calculation at this point.
    if all(k_indexes < (grid_size - n))                 % If we can perform a forward difference calculation... (i.e., if we are not going to go out of bounds.)
        
        % Compute the gradient in each direction at this location.
        gradfk = FowardDifferenceAtPoint(fs, grid, n, k);
        
        % Reshape the gradient in each direction at this location so that it can be stored.
        gradfk = reshape( gradfk, [ones(1, length(num_dims)), num_dims] );
        
        % Create the indexes at which to store this gradient calculation.
        assign_indexes = [num2cell(k_indexes), {':'}];
        
        % Store the gradient at this location into the gradient array.
        gradf = IndexArrayWithArray(gradf, assign_indexes, gradfk);
    end
    
end

end

