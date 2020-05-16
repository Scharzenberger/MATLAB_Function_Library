function k_stimulus = StimulusTensionCurve( r_stimulus, pr, pk, r_offset, k_offset )

% This function computes the stimulus factor k_stimulus used, along with the length factor and maximum active muscle force, to compute the active force in a muscle.

% Inputs:
    % r_stimulus = Scalar Stimulus Ratio.  r_stimulus = Stimulus/Max Stimulus.
    % px = Scalar x convergence ratio.  Determines how close the stimulus is to its maximum value when the tension reaches py of its maximum value.
    % py = Scalar y convergence ratio.  Determines how close the tension is to its maximum value when the stimulus reaches px of its maximum value.
    % r_offset = Scalar r offset.  Determines how much to shift the stimulus tension curve to the right.
    % k_offset = Scalar k offset.  Determines how much to shift the stimulus tension curve up.

% Outputs:
    % k_stimulus = Scalar Stimulus Factor.
    
% Define the default input arguements.
if nargin < 5, k_offset = 0; end
if nargin < 4, r_offset = 0.5; end
if nargin < 3, pk = 0.995; end
if nargin < 2, pr = 0.995; end

% Compute the slope of the sigmoid.
S = (-1./(pr - r_offset)).*log(1./pk - 1);

% Compute the stimulus factor.
k_stimulus = 1./(1 + exp(-S.*(r_stimulus - r_offset))) + k_offset;


end

