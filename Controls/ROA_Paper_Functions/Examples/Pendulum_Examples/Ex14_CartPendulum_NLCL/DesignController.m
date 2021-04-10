%% Design Controller

% This script designs a continuous state space controller with an outer loop intergrator for a simple pendulum linearized about the bottom.

% Clear Everything.
clear, close('all'), clc


%% Define the System Parameters.

% Define the system parameters.
m1_value = 1;
m2_value = 1;
c1_value = 1;
c2_value = 1;
k1_value = 1;
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
% C_bottom = eye(size(A_bottom));
C_bottom = [0 0 1 0];
D_bottom = 0;

% Create the state space model for the cart pendulum linearized at the bottom.
Gs_lbol = ss(A_bottom, B_bottom, C_bottom, D_bottom);

% Create the system matrices for the cart pendulum linearized at the top.
A_top = eval(subs(J_top, [m1 m2 c1 c2 k1 k2 L g F kmotor], [m1_value m2_value c1_value c2_value k1_value k2_value L_value g_value F_value kmotor_value]));
B_top = [0; kmotor_value; 0; 0];
% C_top = eye(size(A_top));
C_top = [0 0 1 0];
D_top = 0;

% Create the state space model for the cart pendulum linearized at the top.
Gs_ltol = ss(A_top, B_top, C_top, D_top);


%% Design a State Space Controller for the Simple Pendulum Linearized About the Bottom.

% Define the desired system response characteristics.
tsettle = 10; PMO = 0;
% tsettle = 5; PMO = 0;
% tsettle = 1; PMO = 0;

% Compute the necessary 2nd order system parameters to achieve the desired system response characteristics.
zeta = PMO2zeta(PMO);
omegan = SettlingTime2omegan(tsettle, zeta);

% Compute the required leading root to achieve the desired 2nd system parameters.
ps = GetsDesignPoints( omegan, zeta );

% Compute the necessary closed loop system roots.
% ps = [ps 10*ps 11*ps 12*ps];
ps = [ps 10*ps 11*ps 12*ps 13*ps];

% Design a state space controller with an outer loop integrator to achieve the desired roots.
% [ Gs_lbcl, K_lbcl ] = GetSSController( Gs_lbol, ps );
[ Gs_lbcl, K_lbcl ] = GetSSIController( Gs_lbol, ps );

% Design a state space controller with an outer loop integrator to achieve the desired roots.
% [ Gs_ltcl, K_ltcl ] = GetSSController( Gs_ltol, ps );
[ Gs_ltcl, K_ltcl ] = GetSSIController( Gs_ltol, ps );


%% Compute the Dynamic Response of the Linearized Cart Pendulum Systems.

% Define the simulation properties for the bottom linearized cart pendulum.
x0_lbol = [0; 0; (pi/180)*45; 0]; x0_lbcl = [0; 0; (pi/180)*45; 0; 0];
tfinal_lbol = 5; tfinal_lbcl = 5;

% Simulate the bottom linearized pendulum.
[ys_lbol, ts_lbol, xs_lbol] = initial(Gs_lbol, x0_lbol, tfinal_lbol);
[ys_lbcl, ts_lbcl, xs_lbcl] = initial(Gs_lbcl, x0_lbcl, tfinal_lbcl);

% Define the simulation properties for the top linearized cart pendulum.
x0_ltol = [0; 0; (pi/180)*45; 0]; x0_ltcl = [0; 0; (pi/180)*45; 0; 0];
tfinal_ltol = 5; tfinal_ltcl = 5;

% Simulate the top linearized pendulum.
[ys_ltol, ts_ltol, xs_ltol] = initial(Gs_ltol, x0_ltol, tfinal_ltol);
[ys_ltcl, ts_ltcl, xs_ltcl] = initial(Gs_ltcl, x0_ltcl, tfinal_ltcl);


%% Plot the Response of the Linearized Cart Pendulums.

% Plot the bottom linearized pendulum results.
figure('Name', 'Bottom Linearized Cart Pendulum (Open Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_lbol, xs_lbol(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_lbol, (180/pi)*xs_lbol(:, 3), 'Linewidth', 3)

figure('Name', 'Bottom Linearized Cart Pendulum (Closed Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_lbcl, xs_lbcl(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_lbcl, (180/pi)*xs_lbcl(:, 3), 'Linewidth', 3)

% Plot the top linearized pendulum  results.
figure('Name', 'Top Linearized Cart Pendulum (Open Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_ltol, xs_ltol(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_ltol, (180/pi)*(xs_ltol(:, 3) + pi), 'Linewidth', 3)

figure('Name', 'Top Linearized Cart Pendulum (Closed Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_ltcl, xs_ltcl(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_ltcl, (180/pi)*(xs_ltcl(:, 3) + pi), 'Linewidth', 3)


