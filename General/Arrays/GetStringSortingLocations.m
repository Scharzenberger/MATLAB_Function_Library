function locs = GetStringSortingLocations(original_cell, desired_cell)

num_original_cells = length(original_cell);
num_desired_cells = length(desired_cell);

indexes = zeros(num_desired_cells, num_original_cells);

for k1 = 1:num_original_cells
    for k2 = 1:num_desired_cells
        
        indexes(k2, k1) = strcmp( original_cell(k1), desired_cell(k2) );
        
    end
end

[locs, ~] = find(indexes');


end

