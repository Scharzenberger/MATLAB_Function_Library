function Keq = GetCableEqStiffness(E, A, L)

% This function computes the equivalent stiffness of a cantilevered beam for displacement at the end of the beam.

% Assumptions:
    % Constant Cross-Section.
    
% Inputs:
    % E = Young's Modulus of Elasticity
    % A = Cross-Sectional Area.
    % L = Cable Length.
    
% Outputs:
    % Keq = Equivalent Spring Stiffness.

% Compute the equivalent spring stiffness.
Keq = A.*E./L;
    
end

