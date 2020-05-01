function Js = T2MOI(Mlist, Ilist, mlist, Ts)

% This function computes the moments of inertia of each joint in an open kinematic chain about in their local frame when each joint is in the position described by Ts (the transformation matrix with respect to the space frame).

% Compute the relative transformation matrices.
Ts_rel = TSpace2TRelative(Ts);

% Compute the moments of inertia of the links in each orientation.
Js = RelT2MOIs(Mlist, Ilist, mlist, Ts_rel);

end

