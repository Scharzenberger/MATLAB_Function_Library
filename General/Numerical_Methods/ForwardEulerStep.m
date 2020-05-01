function U = ForwardEulerStep(U, dU, dt)

% This function performs a single forward Euler step.

% Compute the membrane voltage at the next time step.
U = U + dt*dU;

end

