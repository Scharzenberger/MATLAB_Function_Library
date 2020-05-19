function B = IndexArrayOfUnknownSize(A, dim, index)

% Set the subscripting type to work with normal arrays.
S.type = '()';

% Initialize the subscripting property to select all entries.
S.subs = repmat( {':'}, [1, ndims(A)] );

% Set the subscripting property to select the indexes of interest along the dimension of interest.
S.subs{dim} = index;

% Retrieve the desired part of the array.
B = subsref(A, S);


end

