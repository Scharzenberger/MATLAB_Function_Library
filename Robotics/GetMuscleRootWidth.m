function Lwidth0 = GetMuscleRootWidth(Lwidth, k_width)

% This function computes a muscle's root width, Lwidth0, given the muscle width, Lwidth, and associated length-tension factor, k_width.

% Compute the muscle's root width.
Lwidth0 = Lwidth./sqrt(1 - k_width);


end

