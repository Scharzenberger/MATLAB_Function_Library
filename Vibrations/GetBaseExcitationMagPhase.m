function [mag, phase] = GetBaseExcitationMagPhase(wf, wn, zeta)

% This function computes the magnitude and phase associated with a base excitation mass-spring-damper system.

% Inputs:
    % wf = Forcing Angular Frequency
    % wn = Natural Frequency
    % zeta = Damping Ratio
    
% Outputs:
    % mag = Output/Input Magnitude
    % phase = Output Phase Shift

% Compute the output / input magnitude.
mag = sqrt( (4*(zeta.^2).*((wf./wn).^2) + 1)./( (1 - (wf./wn).^2).^2 + 4*(zeta.^2).*((wf./wn).^2)) );

% Compute the output / input phase.
phase = -atan2( -2*zeta.*((wf./wn).^3), 1 - (wf./wn).^2 + 4*(zeta.^2).*((wf./wn).^2));

end

