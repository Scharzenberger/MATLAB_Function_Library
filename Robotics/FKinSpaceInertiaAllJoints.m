function [Ts, Js] = FKinSpaceInertiaAllJoints(Mlist, Slist, Ilist, mlist, thetalist)

% Retrieve the transformation matrices associated with the given angles.
Ts = FKinSpaceAllJoints(Mlist, Slist, thetalist);

% Compute the moment of inertia of each link in its local frame in each of the robot orientations.
Js = T2MOI(Mlist, Ilist, mlist, Ts);

end

