function J = GetPolarMomentOfHollowRoundBar(d1, d2)

% This function computes the polar second moment of the area of a hollow round bar given its inner and outer diameters d1 and d2.

% Inputs:
    % d1 = Inner Diameter.
    % d2 = Outer Diameter.

% Outputs:
    % J = Polar Second Moment of the Area.
    
% Compute the polar second moment of the area.
J = (pi/32)*(d2.^4 - d1.^4);

end

