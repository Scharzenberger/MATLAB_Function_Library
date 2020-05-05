function epsilon = Length2Strain1D(len, len0)

% This function computes the strain of a member given its current length and starting length.

% Inputs:
    % len = Current length of the member.
    % len0 = Resting length of the member.
    
% Outputs:
    % epsilon = Strain of the member.

% Compute the strain in the member.
epsilon = (len - len0)./len0;

end

