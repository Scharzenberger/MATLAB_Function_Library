function G = Ev2G(E, v)

% This function computes the shear modulus of elasticity given the normal modulus of elasticity and poisson's ratio.

% Assumptions:
    % Linear, Isotropic, Homogeneous Material

% Inputs:
    % E = Young's (Normal) Modulus of Elasticity.  
    % v = Poisson's Ratio.

% Outputs:
    % G = Shear Modulus of Elasticity.

% Compute the shear modulus.
G = E./(2*(1 + v));    

end

