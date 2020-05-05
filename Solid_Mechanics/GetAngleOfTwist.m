function theta = GetAngleOfTwist(T, L, G, J)

% This function computes the angle of twist in a bar due to an applied torsional moment.

% Inputs:
    % T = Applied Torsional Moment.
    % L = Length of Bar.
    % G = Modulus of Rigidity.
    % J = Polar Second Moment of the Area.
    
% Outputs:
    % theta = Angle of Twist.
    
% Compute the angle of twist.
theta = (T.*L)./(G*J);

end

