function [xs, ys, zs, us, vs, ws] = GetQuiverDataFromT(T)

% Retrieve the rotational and translational components associated with this transformation matrix.
[R, P] = TransToRp(T);

% Retrieve the location quiver data.
xs = P(1)*ones(3, 1);
ys = P(2)*ones(3, 1);
zs = P(3)*ones(3, 1);

% Retrieve the direction quiver data.
us = R(1, :)';
vs = R(2, :)';
ws = R(3, :)';

end

