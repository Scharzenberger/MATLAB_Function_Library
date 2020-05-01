function Ls_muscles = GetMuscleLengths(Ps_muscles)

% This function computes the muscle lengths associated with the muscles whose points are described by Ps_muscles

% Create a matrix to store the muscle lengths
Lmuscles_desired = zeros(num_muscles, num_timesteps);

% Compute the muscle lengths throughout the desired trajectory.
for k1 = 1:num_timesteps            % Iterate through each of the time steps...
    for k2 = 1:num_muscles          % Iterate through each of the muscles...
        
        % Compute the distance between the muscle attachment points for this muscle at this time step.
        dPmuscles_desired = diff(Ps_muscles(:, :, k2, k1), 1, 2);
        
        % Compute the length of this muscle at this time step.
        Lmuscles_desired(k2, k1) = sum(vecnorm(dPmuscles_desired, 2, 1));
        
    end
end



end

