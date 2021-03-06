function [ ks, V_offset ] = GetBipolar2SingleEndedGains( V_actual, V_target, V_supply )
%This function computes the gains for a bipolar to single ended op-amp circuit.  The first gain is the gain required to scale the current signal into the desired signal domain.  The second gain is the gain required to shift the current signal into the desired signal domain.

%INPUTS:
    %V_actual = Current voltage domain as a 1x2 array.  eg, [V_actual_low V_actual_high].
    %V_target = Target voltage domain as a 1x2 array.  eg, [V_target_low V_target_high].
    %V_supply = Available supply voltage for the bipolar to single ended offset value.
    
%OUTPUTS:
    %ks = Required op-amp gains as a 1x2 array.
    %V_offset = The voltage offset that is generated by the second required gain.
    
%Compute the gain for the bipolar to single ended op-amp circuit.
k1 = diff(V_target)/diff(V_actual);

%Compute the required voltage offset value for the bipolar to single ended op-amp circuit.
V_offset = k1*V_actual(1);

%Compute the gain necessary to achieve the required voltage offset.
k2 = V_offset/V_supply;

%Store these gains into an array.
ks = [k1 k2];

end

