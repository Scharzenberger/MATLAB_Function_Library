function L = InverseLengthTensionCurve( k_length, Lrest, Lwidth, k_width )

% This function computes the muscle length associated with a specific muscle length factor.

% Inputs:
    % k_length = Scalar Length Factor.
    % Lrest = Scalar Resting Muscle Length.
    % Lwidth = Scalar Muscle "Width" at which the Length Factor, k_length, becomes k_width. i.e., when L = Lrest +- Lwidth, then k_length = k_width.
    % k_width = Scalar Length Factor that Occurs When L = Lrest +- Lwidth.

% Outputs:
    % L = 2x1 Vector of Muscle Lengths.

% Set the default input arguments.
if nargin < 4, k_width = 0; end       % Note that if k_width = 0, then Lwidth0 = Lwidth.  i.e., Lwidth is the muscle "width" at which the muscle produces zero force (since k_length = 0).
    
% Compute the width of the parabola to its roots.
Lwidth0 = Lwidth./sqrt(1 - k_width);

% Compute the muscle length.
L = [Lrest - Lwidth0.*sqrt(1 - k_length); Lrest + Lwidth0.*sqrt(1 - k_length)];

end

