function r_stimulus = InverseStimulusTensionCurve( k_stimulus, pr, pk, r_offset, k_offset )

% This function computes the stimulus ratio associated with a specific stimulus factor.

% Inputs:
    % k_stimulus = Scalar Stimulus Factor.
    % pr = Scalar r convergence ratio.  Determines how close the stimulus is to its maximum value when the tension reaches pk of its maximum value.
    % pk = Scalar k convergence ratio.  Determines how close the tension is to its maximum value when the stimulus reaches pr of its maximum value.
    % r_offset = Scalar r offset.  Determines how much to shift the stimulus tension curve to the right.
    % k_offset = Scalar k offset.  Determines how much to shift the stimulus tension curve up.

% Outputs:
    % r_stimulus = Scalar Stimulus Ratio.  stim_ratio = Stimulus/Max Stimulus.

% Define the default input arguements.
if nargin < 5, k_offset = 0; end
if nargin < 4, r_offset = 0.5; end
if nargin < 3, pk = 0.995; end
if nargin < 2, pr = 0.995; end

% Compute the slope of the sigmoid.
S = (-1./(pr - r_offset)).*log(1./pk - 1);

% Compute the stimulus ratio.
r_stimulus = -(1./S).*log( 1./(k_stimulus - k_offset) - 1 ) + r_offset;


end

