%% Design Controller

% This script designs a continuous state space controller with an outer loop intergrator for a simple pendulum linearized about the bottom.

% Clear Everything.
clear, close('all'), clc


%% Define the Simple Pendulum Linearized About The Bottom.

% Define the pendulum parameters.
m = 1; L = 1; g = 9.81; c = 1; ktau = 0; kmotor = 1;

% Define the linearized system matrices.
A = [0 1; -(ktau + m*g*L)/(m*(L^2)) -c/(m*(L^2))];
B = [0; kmotor/(m*(L^2))];
C = [1 0];
D = [0];

% Define the open loop simple pendulum linearized about the bottom.
Gol = ss(A, B, C, D);


%% Design a State Space Controller for the Simple Pendulum Linearized About the Bottom.

% Define the desired system response characteristics.
tsettle = 1; PMO = 0;

% Compute the necessary 2nd order system parameters to achieve the desired system response characteristics.
zeta = PMO2zeta(PMO);
omegan = SettlingTime2omegan(tsettle, zeta);

% Compute the required leading root to achieve the desired 2nd system parameters.
ps = GetsDesignPoints( omegan, zeta );

% Compute the necessary closed loop system roots.
ps = [ps 10*ps 11*ps];

% Design a state space controller with an outer loop integrator to achieve the desired roots.
[ Gcl, K ] = GetSSIController( Gol, ps );


