function B = IndexArrayWithArray(A, indexes, V)

% Set the index type to work for arrays.
S.type = '()';

% Set the indexes to reference to be those specified in indexes.
if iscell(indexes)                      % If the indexes are already stored as a cell array...
    
    % Directly assign the indexes into our substitution variable.
    S.subs = indexes; 
    
else                                    % Otherwise...
    
    % Convert the indexes to a cell array before storing the indexes into our substitution variable.
    S.subs = num2cell(indexes);
    
end

% Determine whether to retrieve the indexes from A or replace them with V.
if nargin < 3                   % If there was no specified value to assign...
    
    % Retrieve the desired entries from A.
    B = subsref(A, S);
    
else                            % Otherwise...
    
    % Store the desired entries into A.
    B = subsasgn(A, S, V);

end

end

