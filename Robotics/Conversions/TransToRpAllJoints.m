function [Rs, Ps] = TransToRpAllJoints(Ts)

% This function retrieves the rotational and translational components associated with the given transformation matrices.

% Retrieve the number of transformation matrices.
num_joints = size(Ts, 3);
num_angles = size(Ts, 4);

% Initialize tensors to store the rotational and translational components that we will extract from the transfomratino matrices.
Rs = zeros(3, 3, num_joints, num_angles);
Ps = zeros(3, num_joints, num_angles);

% Retrieve the rotational and translational components associated with each transformation matrix.
for k1 = 1:num_angles                   % Iterate through each of the angles...
    for k2 = 1:num_joints                % Iterate through each of the joints...
        
        % Retrieve the rotational and translational components associated with this transformation matrix.
        [Rs(:, :, k2, k1), Ps(:, k2, k1)] = TransToRp(Ts(:, :, k2, k1));
        
    end
end

end

