function I = TransformMomentOfInertia(I, m, T)

% This function compute the moment of inertia J after being transformed by T.

% Retrieve the rotationa and translational components of the transformation matrix.
[R, P] = TransToRpAllJoints(T);

% Rotate the inertia matrix.
I = RotateMomentOfInertia(I, R);

% Translate the inertia matrix.
I = TranslateMomentOfInertia(I, m, P);

end

