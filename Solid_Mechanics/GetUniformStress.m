function sigma = GetUniformStress(F, A)

% This function computes the uniform stress caused by a force F over an area A.

% Inputs:
    % F = Uniform Force
    % A = Cross-Sectional Area

% Outputs:
    % sigma = Uniform Stress
    
% Compute the uniform stress.
sigma = F./A;


end

