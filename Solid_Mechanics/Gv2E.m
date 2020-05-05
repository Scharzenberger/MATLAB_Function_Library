function E = Gv2E(G, v)

% This function computes the normal modulus of elasticity given the shear modulus of elasticity and poisson's ratio.

% Assumptions:
    % Linear, Isotropic, Homogeneous Material

% Inputs:
    % G = Shear Modulus of Elasticity.
    % v = Poisson's Ratio.

% Outputs:
    % E = Young's (Normal) Modulus of Elasticity.

% Compute Young's modulus.
E = 2*G*(1 + v);
    

end

