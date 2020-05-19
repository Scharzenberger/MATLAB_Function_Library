function gradf = FowardDifferenceAtPointInDirection(fs, grid, n, loc, dir)

% Collect size information.
dims = size(fs);

% Retrieve the grid points associated with this direction.
grid_slice = IndexArrayOfUnknownSize(grid, ndims(grid), dir);

% Define the starting index values.
grid_slice_indexes0 = AbsInd2DimInd( dims, loc );
f_indexes0 = AbsInd2DimInd( dims, loc ); f_indexes0(dir) = f_indexes0(dir) + n;

% Initialize the starting index values.
grid_slice_indexes_low = grid_slice_indexes0; grid_slice_indexes_high = grid_slice_indexes0;
f_indexes = f_indexes0;

% Initialize the dx and df calculations.
dx = 1;
df = IndexArrayWithArray(fs, f_indexes0);

% Compute the dx and df approximations via an nth order foward difference method.
for i = 1:n             % Iterate through each of the orders...
    
    % Advance the x index variable.
    grid_slice_indexes_low(dir) = grid_slice_indexes0(dir) + i - 1;
    grid_slice_indexes_high(dir) = grid_slice_indexes0(dir) + i;

    % Advance the f index variable.
    f_indexes(dir) = f_indexes0(dir) - i;
    
    % Retrieve the low and high x values.
    x_low = IndexArrayWithArray(grid_slice, grid_slice_indexes_low);
    x_high = IndexArrayWithArray(grid_slice, grid_slice_indexes_high);
    
    % Approximate the dx term.
    dx = dx.*(x_high - x_low);

    % Retrieve the relevant function value.
    fk = IndexArrayWithArray(fs, f_indexes);
    
    % Approximate the df term.
    df = df + ((-1).^i).*nchoosek(n, i).*fk;

end

% Approximate the nth order df/dx calculation.
gradf = df./dx;


end

