function gradf = FowardDifferenceAtPoint(fs, grid, n, loc)

% Retrieve the number of dimensions from the input sizes.
num_dims = ndims(fs);

% Preallocate an array to store the gradient at this location.
gradf = zeros(num_dims, 1);

% Estimate the derivative in each direction at the given location.
for k = 1:num_dims                  % Iterate through each dimension (i.e., direction)...

    % Estimate the derivative at this point in this direction.
    gradf(k) = FowardDifferenceAtPointInDirection(fs, grid, n, loc, k);

end

end

