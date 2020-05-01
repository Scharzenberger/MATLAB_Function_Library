function Ts = FKinSpaceAllJoints(Mlist, Slist, thetalist)

% This functions computes the transformation matrix associated with each joint in an open kinematic chain.

% Retrieve information about the size of our input arguments.
num_joints = size(thetalist, 1);
num_angles = size(thetalist, 2);
num_Ss = size(Slist, 2);
num_Ms = size(Mlist, 3);

% Ensure that there are the same number of screw axes as given joint angles.
if (num_Ss ~= num_joints)
    
    % Throw an error stating that the number of screw axes must agree with the number of provided joint angles.
    error('Number of screw axes (size(Slist, 2) = %0.0f) must equal the number of provided joint angles (length(thetalist) = %0.0f).', size(Slist, 2), num_joints)
    
end

% Ensure that there are the same number of home transformation matrices as provided joint angles (or one more transformation matrix than joint angles).
if (num_Ms ~= num_joints) && (num_Ms ~= (num_joints + 1))
    
    % Throw an error stating that the number
    error('Number of home matrices (size(Mlist, 3) = %0.0f) must equal either: \n(1) the number of provided joint angles (length(thetalist) = %0.0f), OR \n(2) the number of provided joint angles plus one (length(thetalist) + 1 = %0.0f)', num_Ms, num_joints, num_joints + 1)
    
end

% Determine whether we have a grounded point.
if num_Ms > num_joints          % If the number of home matrices is greater than the number of joints...
    
    % Retrieve the stationary home matrix.
    M0 = Mlist(:, :, 1);
    
    % Set the home matrix list equal to all but the first home matrix.
    Mlist = Mlist(:, :, 2:end);
    
else
    
    % Set the stationary home matrix to be empty.
    M0 = [];
    
end

% Initialize a matrix to store the transformation matrix associated with each of the joints.
Ts = zeros(4, 4, num_joints + 1, num_angles);

% Compute the transformation matrix associated with each angle and each joint.
for k1 = 1:num_angles                       % Iterate through each of the angles...
    
    for k2 = 1:num_joints                   % Iterate through each of the joints...
        
        % Compute the transformation matrix associated with the current joint and the current angles.
        Ts(:, :, k2 + 1, k1) = FKinSpace(Mlist(:, :, k2), Slist(:, 1:k2), thetalist(1:k2, k1));
        
    end
    
    Ts(:, :, 1, k1) = Ts(:, :, 2, k1);
    
    Ts(1:3, 4, 1, k1) = zeros(3, 1);
    
end

end

