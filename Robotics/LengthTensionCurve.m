function [k_length, Lwidth0] = LengthTensionCurve( L, Lrest, Lwidth, k_width )

% This function computes the length factor used, along with the stimulus factor and maximum active muscle force, to compute the active force in a muscle.

% Inputs:
    % L = Scalar Muscle Length.
    % Lrest = Scalar Resting Muscle Length.
    % Lwidth = Scalar Muscle "Width" at which the Length Factor, k_length, becomes k_width. i.e., when L = Lrest +- Lwidth, then k_length = k_width.
    % k_width = Scalar Length Factor that Occurs When L = Lrest +- Lwidth.

% Outputs:
    % k_length = Scalar Length Factor.
    
% Set the default input arguments.
if nargin < 4, k_width = 0; end       % Note that if k_width = 0, then Lwidth0 = Lwidth.  i.e., Lwidth is the muscle "width" at which the muscle produces zero force (since k_length = 0).
    
% Compute the width of the parabola to its roots.
Lwidth0 = GetMuscleRootWidth(Lwidth, k_width);

% Compute the length factor.
k_length = 1 - ((L - Lrest).^2)./(Lwidth0.^2);

end

