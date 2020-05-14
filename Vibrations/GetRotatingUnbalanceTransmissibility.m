function TR = GetRotatingUnbalanceTransmissibility(wf, wn, zeta)

% This function computes the transmissibility associated with a rotating unbalance mass-spring-damper system.

% Inputs:
    % wf = Forcing Angular Frequency
    % wn = Natural Frequency
    % zeta = Damping Ratio
    
% Outputs:
    % TR = Transmissibility

% Compute the transmissibility.
TR = sqrt( (4*(zeta.^2).*((wf./wn).^2) + 1)./( (1 - (wf./wn).^2).^2 + 4*(zeta.^2).*((wf./wn).^2) ) );


end

