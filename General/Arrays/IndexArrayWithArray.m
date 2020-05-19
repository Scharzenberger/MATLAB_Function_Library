function B = IndexArrayWithArray(A, indexes)

% Set the index type to work for arrays.
S.type = '()';

% Set the indexes to reference to be those specified in indexes.
S.subs = num2cell(indexes);

% Retrieve the desired entries from A.
B = subsref(A, S);


end

