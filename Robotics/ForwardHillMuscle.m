function [T, dT] = ForwardHillMuscle(T0, L, dL, A, kse, kpe, b, dt, intRes)

% NOT CURRENTLY USING INTRES, BUT POTENTIALLY NEED TO.

% This function computes the muscle tension and rate of change of the muscle tension over time in a Hill Muscle given an activation and length history.

% Inputs:
% T0 = Initial Muscle Tension as a num_muscles x 1 vector.
% L = Muscle Length as a num_muscles x num_timesteps matrix.
% dL = Muscle Velocity as a num_muscles x num_timesteps matrix.
% A = Muscle Activation as a num_muscles x num_timesteps matrix.
% kse = Series Muscle Stiffness.
% kpe = Parallel Muscle Stiffness.
% b = Damping Coefficient.

% Outputs:
% T = Muscle Tension as a num_muscles x num_timesteps matrix.
% dT = Rate of Change of Muscle Tension Over Time as a num_muscles x num_timesteps matrix.

% Retrieve size information from the inputs.
num_muscles = size(L, 1);
num_timesteps = size(L, 2);

% Initialize matrices to store the muscle tensions.
[T, dT] = deal( zeros(num_muscles, num_timesteps) );

% Set the initial tension.
T(:, 1) = T0;

% Compute the muscle tensions & rate of change of muscle tension over time.
for k = 1:(num_timesteps - 1)

    % Compute the rate of change of the muscle tension at this time step.
    dT(:, k) = (kse/b)*(kpe*L(:, k) + b*dL(:, k) - (1 + (kpe/kse))*T(:, k) + A(:, k));

    % Compute the muscle tension at the next time step.
    T(:, k + 1) = T(:, k) + dt*dT(:, k);

end

% Compute the final rate of change of the muscle tension.
dT(:, k + 1) = (kse/b)*(kpe*L(:, k + 1) + b*dL(:, k + 1) - (1 + (kpe/kse))*T0 + A(:, k + 1));


% % Compute the integration step size.
% dt = dt/intRes;
% 
% % Compute the muscle tensions & rate of change of muscle tension over time.
% for k1 = 1:(num_timesteps - 1)
%     
%     dTint = (kse/b)*(kpe*L(:, k1) + b*dL(:, k1) - (1 + (kpe/kse))*T0 + A(:, k1));
%     
%     for k2 = 1:intRes
%         
%         % Compute the rate of change of the muscle tension at this time step.
%         dT(:, k1) = (kse/b)*(kpe*L(:, k1) + b*dL(:, k1) - (1 + (kpe/kse))*T0 + A(:, k1));
%         
%         % Compute the muscle tension at the next time step.
%         Tint = Tint + dt*dTint;
%         
%     end
%     
%     
%     % Compute the rate of change of the muscle tension at this time step.
%     dT(:, k1) = dTint;
%     
%     % Compute the muscle tension at the next time step.
%     T(:, k1 + 1) = T(:, k1) + dt*dT(:, k1);
%     
%     
% end
% 
% % Compute the final rate of change of the muscle tension.
% dT(:, k1 + 1) = (kse/b)*(kpe*L(:, k1 + 1) + b*dL(:, k1 + 1) - (1 + (kpe/kse))*T0 + A(:, k1 + 1));


end

