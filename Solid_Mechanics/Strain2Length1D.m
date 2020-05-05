function len = Strain2Length1D(epsilon, len0)

% This function computes the strain of a member given its current length and starting length.

% Inputs:
    % epsilon = Strain of the member.
    % len0 = Resting length of the member.
    
% Outputs:
    % len = Current length of the member.

% Compute the strain in the member.
len = (1 + epsilon).*len0;

end

