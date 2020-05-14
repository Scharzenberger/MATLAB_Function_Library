function v = EG2v(E, G)

% This function computes poisson's ratio given the normal and shear modulii of elasticity.

% Assumptions:
    % Linear, Isotropic, Homogeneous Material

% Inputs:
    % E = Young's (Normal) Modulus of Elasticity.  
    % G = Shear Modulus of Elasticity.

% Outputs:
    % v = Poisson's Ratio.

% Compute Poisson's ratio.
v = E./(2*G) - 1;

end

