function J = GetPolarMomentOfSolidRoundBar(d)

% This function computes the polar second moment of the area of a solid round bar given its diameter d.

% Inputs:
    % d = Diameter of the Solid Round Bar.

% Outputs:
    % J = Polar Second Moment of the Area.
    
% Compute the polar second moment of the area.
J = (pi/32)*(d.^4);


end

