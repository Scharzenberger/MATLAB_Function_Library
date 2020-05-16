function [Factive, k_stimulus, k_length] = LengthStimulus2ActiveForce( Factive_max, r_stimulus, L, Lrest, Lwidth, k_width, pr, pk, r_offset, k_offset )

% This function computes the active force in a muscle given its current length, stimulus ratio, and maximum active force.

% Inputs:
    % Factive_max = Scalar Maximum Active Muscles Force.
    % r_stimulus = Scalar Stimulus Ratio.  r_stim = Stimulus/Max Stimulus.
    % Lrest = Scalar Resting Muscle Length.
    % Lwidth = Scalar Muscle "Width" at which the Length Factor, k_length, becomes k_width. i.e., when L = Lrest +- Lwidth, then k_length = k_width.
    % k_width = Scalar Length Factor that Occurs When L = Lrest +- Lwidth.
    % px = Scalar x convergence ratio.  Determines how close the stimulus is to its maximum value when the tension reaches py of its maximum value.
    % py = Scalar y convergence ratio.  Determines how close the tension is to its maximum value when the stimulus reaches px of its maximum value.
    % r_offset = Scalar r offset.  Determines how much to shift the stimulus tension curve to the right.
    % k_offset = Scalar k offset.  Determines how much to shift the stimulus tension curve up.

% Outputs:
    % Factive = Scalar Active Muscle Force.
    % k_stimulus = Scalar Stimulus Factor.
    % k_length = Scalar Length Factor.

% Define the default input arguements.
if nargin < 10, k_offset = 0; end
if nargin < 9, r_offset = 0.5; end
if nargin < 8, pk = 0.995; end
if nargin < 7, pr = 0.995; end
if nargin < 6, k_width = 0; end       % Note that if k_width = 0, then Lwidth0 = Lwidth.  i.e., Lwidth is the muscle "width" at which the muscle produces zero force (since k_length = 0).

% Compute the length factor.
k_length = LengthTensionCurve( L, Lrest, Lwidth, k_width );
    
% Compute the stimulus factor.
k_stimulus = StimulusTensionCurve( r_stimulus, pr, pk, r_offset, k_offset );
    
% Compute the active force.
Factive = k_stimulus.*k_length.*Factive_max;


end

