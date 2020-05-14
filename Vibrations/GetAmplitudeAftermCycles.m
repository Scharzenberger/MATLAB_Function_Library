function Am = GetAmplitudeAftermCycles(A0, zeta, m)

% This function computes the amplitude ym of a mass-spring-damper system after m cycles given the starting amplitude y0 and damping ratio zeta.

% Inputs:
    % A0 = Initial Amplitude
    % zeta = Damping Ratio
    % m = Number of Cycles.
    
% Outputs:
    % Am = Amplitude After m Cycles.

% Compute the amplitude after m cycles.
Am = A0.*exp( -(2*pi*m.*zeta)./sqrt(1 - zeta.^2) );

end

