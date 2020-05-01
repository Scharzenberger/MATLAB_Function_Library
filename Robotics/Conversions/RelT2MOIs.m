function Jlist = RelT2MOIs(Mlist, Ilist, mlist, Ts_rel)

% This function computes the moments of inertia associated with the links of an open kinematic chain Jlist at their local joints in their local frame for the given relative transformation matrices given their moments of interia about their local joints in their local frame in the home position.

% Inputs:
% Ilist = High dimensional array that includes the inertia matrices of each joint about their local joint in their local frame when the open kibmeatic chain is in the home orientation.
% Ts_rel = High dimensiondal array that includes the relative transformation matrices for each joint relative to their most proximal neighboring joint.

% Outputs:
% Jlist = High dimensional array that includes the inertia matrices of each joint about their local joint in their local frame when the open kinematic chain is in the orientation described by Ts_rel.

% Retrieve size information from the inputs.
num_joints = size(Ts_rel, 3);
num_angles = size(Ts_rel, 4);

% Conver the Mlist into a relative Mlist.
Mlist_rel = TSpace2TRelative(Mlist);

% Preallocate a matrix to store the new moments of inertia.
Jlist = zeros([size(Ilist) num_angles]);

% Compute the moment of inertia for each joint in the given orientation.
for k1 = 1:num_angles                       % Iterate through each of the angles...
    
    % Transfer over the moment of inertia associated with the final link (it is not affected by the rotation of the joints).
    Jlist(:, :, end, k1) = Ilist(:, :, end);
    
    for k2 = (num_joints - 1):-1:1             % Iterate through each link less one (i.e., each joint less one)...
        
        % Retrieve the rotational component of the inertia transformation matrix.
        R = Ts_rel(1:3, 1:3, k2 + 1, k1);
        
        % Compute the translational component of the inertia transformation matrix.
        P = -Mlist_rel(1:3, 4, k2);
        
        % Compute the transformation matrix necessary to transform the inertia matrix.
%         T = [Ts_rel(1:3, 1:3, k2 + 1, k1), -Ts_rel(1:3, 4, k2, k1); zeros(1, 3), 1];
        T = [R, P; zeros(1, 3), 1];

        % Compute the moment of inertia of the current link.
        Jlist(:, :, k2, k1) = Ilist(:, :, k2) + TransformMomentOfInertia(Jlist(:, :, k2 + 1, k1), mlist(k2 + 1), T);
        
    end
    
end

end

