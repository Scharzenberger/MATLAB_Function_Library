%% Get Linearized System

% This script computes a linearized system of equations for a cart pendulum system.

% Clear Everything.
clear, close('all'), clc

% Set the default figure color.
set(0,'DefaultFigureColor', 'w');


%% Define the System Parameters.

% Define the system parameters.
m1_value = 1;
m2_value = 1;
c1_value = 1;
c2_value = 1;
k1_value = 0;
k2_value = 0;
L_value = 1;
g_value = 9.81;
F_value = 0;
kmotor_value = 1;


%% Define the Cart Pendulum System of Equations

% Define the symbolic variables.
syms m1 m2 c1 c2 k1 k2 L g F kmotor x1 x2 x3 x4 x1dot x2dot x3dot x4dot

% Define the system of equations.
eq1 = (m1 + m2)*x2dot + c1*x2 + k1*x1 + (1/2)*m2*L*x4*sin(x3) - (1/2)*m2*L*x4dot*cos(x3) == F;
eq2 = ((7/12)*m2*(L^2))*x4dot + c2*x4 + k2*x3 - (1/2)*m2*g*L*sin(x3) - (1/2)*m2*L*cos(x3)*x2dot == 0;

% eq1 = (m1 + m2)*x2dot + c1*x2 - (1/2)*m2*L*x4dot == F;
% eq2 = ((7/12)*m2*(L^2))*x4dot + c2*x4 - (1/2)*m2*g*L*x3 - (1/2)*m2*L*x2dot == 0;

% Solve the system of equations.
sol = solve([eq1, eq2], [x2dot, x4dot]);

% Retrieve the solved system of equations.
x1dot = x2;
x2dot = sol.x2dot;
x3dot = x4;
x4dot = sol.x4dot;

% Define the ODE system.
f = [x1dot; x2dot; x3dot; x4dot];


%% Compute the System Jacobian.

% Set the equilibrium point.
xeq_bottom = [0; 0; pi; 0];
xeq_top = [0; 0; 0; 0];

% Compute the system jacobian.
J = [diff(x1dot, x1) diff(x1dot, x2) diff(x1dot, x3) diff(x1dot, x4); diff(x2dot, x1) diff(x2dot, x2) diff(x2dot, x3) diff(x2dot, x4); diff(x3dot, x1) diff(x3dot, x2) diff(x3dot, x3) diff(x3dot, x4); diff(x4dot, x1) diff(x4dot, x2) diff(x4dot, x3) diff(x4dot, x4)];

% Evaluate the jacobian at the equilibrium point.
J_bottom = subs(J, [x1, x2, x3, x4], xeq_bottom');

% Evaluate the jacobian at the equilibrium point.
J_top = subs(J, [x1, x2, x3, x4], xeq_top');




%% Create a State Space Models of the Linearized Cart Pendulum.

% Create the system matrices for the cart pendulum linearized at the bottom.
A_bottom = eval(subs(J_bottom, [m1 m2 c1 c2 k1 k2 L g F kmotor], [m1_value m2_value c1_value c2_value k1_value k2_value L_value g_value F_value kmotor_value]));
B_bottom = [0; kmotor_value; 0; 0];
C_bottom = [1 0 0 0; 0 0 1 0];
D_bottom = 0;

% Create the state space model for the cart pendulum linearized at the bottom.
Gs_lbol = ss(A_bottom, B_bottom, C_bottom, D_bottom);

% Create the system matrices for the cart pendulum linearized at the top.
A_top = eval(subs(J_top, [m1 m2 c1 c2 k1 k2 L g F kmotor], [m1_value m2_value c1_value c2_value k1_value k2_value L_value g_value F_value kmotor_value]));
B_top = [0; kmotor_value; 0; 0];
C_top = [1 0 0 0; 0 0 1 0];
D_top = 0;

% Create the state space model for the cart pendulum linearized at the top.
Gs_ltol = ss(A_top, B_top, C_top, D_top);


%% Compute the Dynamic Response of the Linearized Cart Pendulum Systems.

% Define the simulation properties for the bottom linearized cart pendulum.
x0_lbol = [0; 0; pi; 0];
tfinal_lbol = 10;

% Simulate the bottom linearized pendulum.
[ys_lbol, ts_lbol] = initial(Gs_lbol, x0_lbol, tfinal_lbol);

% Define the simulation properties for the top linearized cart pendulum.
x0_ltol = [0; 0; pi; 0];
tfinal_ltol = 10;

% Simulate the top linearized pendulum.
[ys_ltol, ts_ltol] = initial(Gs_ltol, x0_ltol, tfinal_ltol);


%% Plot the Response of the Linearized Cart Pendulums.

% Plot the bottom linearized pendulum  results.
figure('Name', 'Bottom Linearized Pendulum')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_lbol, ys_lbol(:, 1))
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s'), ylabel('Pendulum Angle [rad]'), title('Pendulum Angle vs Time'), plot(ts_lbol, ys_lbol(:, 2))

% Plot the top linearized pendulum  results.
figure('Name', 'Top Linearized Pendulum')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_ltol, ys_ltol(:, 1))
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s'), ylabel('Pendulum Angle [rad]'), title('Pendulum Angle vs Time'), plot(ts_ltol, ys_ltol(:, 2))

