%% Design Controller

% This script designs a continuous state space controller with an outer loop intergrator for the quansar pendulum.

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
L1_value = 1;
L2_value = 1;
g_value = 9.81;
etag_value = 1;
etam_value = 1;
kt_value = 1;
km_value = 1;
Kg_value = 1;
Rm_value = 1;
Vm_value = 0;


%% Define the Cart Pendulum System of Equations

% Define the symbolic variables.
syms m1 m2 c1 c2 k1 k2 L1 L2 g etag etam kt km Kg Rm Vm x1 x2 x3 x4 x1dot x2dot x3dot x4dot

% Define the system of equations.
eq1 = ( ((1/12)*m1 + m2)*(L1^2) + (1/4)*m2*(L2^2)*(1 - (cos(x3)^2)) )*x2dot + ( -(1/2)*m2*L2*L1*cos(x3) )*x4dot + ( (1/2)*m2*(L2^2)*sin(x3)*cos(x3) )*x2*x4 + ( (1/2)*m2*L2*L1*sin(x3) )*(x4^2) + ( (etag*etam*kt*km*(Kg^2) + c1*Rm)/Rm )*x2 + k1*x1 == (etag*etam*kt*Kg/Rm)*Vm;
eq2 = ( -(1/2)*m2*L2*L1*cos(x3) )*x2dot + ( (1/3)*m2*(L2^2) )*x4dot + ( -(1/4)*m2*L2^2*cos(x3)*sin(x3) )*(x2^2) + ( (1/2)*m2*L2*g*sin(x3) ) + c2*x4 + k2*x3 == 0;

% Solve the system of equations.
sol = solve([eq1, eq2], [x2dot, x4dot]);

% Retrieve the solved system of equations.
x1dot = x2;
x2dot = sol.x2dot; x2dot = simplify(x2dot); % x2dot = collect(x2dot, [x1, x2, x3, x4]);
x3dot = x4;
x4dot = sol.x4dot; x4dot = simplify(x4dot); % x4dot = collect(x4dot, [x1, x2, x3, x4]);

% Define the ODE system.
f = [x1dot; x2dot; x3dot; x4dot];

% % Compute the system derivatives.
% x1dot = x2;
% x2dot = (9*L1.*L2.*m2.*cos(x3).^2.*sin(x3).*x2.^2)./(2*L1.^2.*m1 + 24*L1.^2.*m2 + 6*L2.^2.*m2 - 18*L1.^2.*m2.*cos(x3).^2 - 6*L2.^2.*m2.*cos(x3).^2) - (6*L2.^2.*m2.*cos(x3).*sin(x3).*x2.*x4)./(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2) - ((24*L2.*etag.*etam.*km.*kt.*Kg.^2 + 24*L2.*Rm.*c1).*x2)./(2*L2.*Rm.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2)) - (6*L1.*L2.*m2.*sin(x3).*x4.^2)./(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2) - (18*L1.*c2.*cos(x3).*x4)./(L2.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2)) + (24*Kg.*L2.*Vm.*etag.*etam.*kt + 18*L1.*L2.*Rm.*g.*m2.*cos(x3).*sin(x3))./(2*L2.*Rm.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2));
% x3dot = x4;
% x4dot = ((3*(12*Rm.*sin(x3).*L1.^2.*L2.^2.*m2.^2.*cos(x3) + Rm.*m1.*sin(x3).*L1.^2.*L2.^2.*m2.*cos(x3) - 3*Rm.*sin(x3).*L2.^4.*m2.^2.*cos(x3).^3 + 3*Rm.*sin(x3).*L2.^4.*m2.^2.*cos(x3)))./(4*L2.^2.*Rm.*m2.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2))).*x2.^2 + (-(9*L1.*L2.*m2.*cos(x3).^2.*sin(x3))./(L1.^2.*m1 + 12*L1.^2..*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2)).*x2.*x4 + (-(3*(24*L1.*L2.*etag.*etam.*km.*kt.*m2.*cos(x3).*Kg.^2 + 24*L1.*L2.*Rm.*c1.*m2.*cos(x3)))./(4*L2.^2.*Rm.*m2.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2))).*x2 + (-(9*L1.^2.*m2.*cos(x3).*sin(x3))./(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2)).*x4.^2 + (-(3*(4*L1.^2.*Rm.*c2.*m1 + 48*L1.^2.*Rm.*c2.*m2 + 12*L2.^2.*Rm.*c2.*m2 - 12*L2.^2.*Rm.*c2.*m2.*cos(x3).^2))./(4*L2.^2.*Rm.*m2.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2))).*x4 + (3*(24*Rm.*g.*sin(x3).*L1.^2.*L2.*m2.^2 + 2*Rm.*g.*m1.*sin(x3).*L1.^2.*L2.*m2 + 24*Kg.*Vm.*etag.*etam.*kt.*L1.*L2.*m2.*cos(x3) - 6*Rm.*g.*sin(x3).*L2.^3.*m2.^2.*cos(x3).^2 + 6*Rm.*g.*sin(x3).*L2.^3.*m2.^2))./(4*L2.^2.*Rm.*m2.*(L1.^2.*m1 + 12*L1.^2.*m2 + 3*L2.^2.*m2 - 9*L1.^2.*m2.*cos(x3).^2 - 3*L2.^2.*m2.*cos(x3).^2));
% 
% % Define the ODE system.
% f = [x1dot; x2dot; x3dot; x4dot];


%% Compute the System Jacobian.

% Set the equilibrium point.
xeq_bottom = [0; 0; 0; 0];
xeq_top = [0; 0; pi; 0];

% Compute the system jacobian.
J = [diff(x1dot, x1) diff(x1dot, x2) diff(x1dot, x3) diff(x1dot, x4); diff(x2dot, x1) diff(x2dot, x2) diff(x2dot, x3) diff(x2dot, x4); diff(x3dot, x1) diff(x3dot, x2) diff(x3dot, x3) diff(x3dot, x4); diff(x4dot, x1) diff(x4dot, x2) diff(x4dot, x3) diff(x4dot, x4)];

% Evaluate the jacobian at the equilibrium point.
J_bottom = subs(J, [x1, x2, x3, x4], xeq_bottom');

% Evaluate the jacobian at the equilibrium point.
J_top = subs(J, [x1, x2, x3, x4], xeq_top');


%% Create a State Space Models of the Linearized Cart Pendulum.

% Create the system matrices for the cart pendulum linearized at the bottom.
A_bottom = eval(subs(J_bottom, [m1 m2 c1 c2 k1 k2 L1 L2 g etag etam kt km Kg Rm Vm], [m1_value m2_value c1_value c2_value k1_value k2_value L1_value L2_value g_value etag_value etam_value kt_value km_value Kg_value Rm_value Vm_value]));
B_bottom = [0; 1; 0; 0];
C_bottom = [0 0 1 0];
D_bottom = 0;

% Create the state space model for the cart pendulum linearized at the bottom.
Gs_lbol = ss(A_bottom, B_bottom, C_bottom, D_bottom);

% Create the system matrices for the cart pendulum linearized at the top.
A_top = eval(subs(J_top, [m1 m2 c1 c2 k1 k2 L1 L2 g etag etam kt km Kg Rm Vm], [m1_value m2_value c1_value c2_value k1_value k2_value L1_value L2_value g_value etag_value etam_value kt_value km_value Kg_value Rm_value Vm_value]));
B_top = [0; 1; 0; 0];
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
figure('Color', 'w', 'Name', 'Bottom Linearized Cart Pendulum (Open Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_lbol, xs_lbol(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_lbol, (180/pi)*xs_lbol(:, 3), 'Linewidth', 3)

figure('Color', 'w', 'Name', 'Bottom Linearized Cart Pendulum (Closed Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_lbcl, xs_lbcl(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_lbcl, (180/pi)*xs_lbcl(:, 3), 'Linewidth', 3)

% Plot the top linearized pendulum  results.
figure('Color', 'w', 'Name', 'Top Linearized Cart Pendulum (Open Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_ltol, xs_ltol(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_ltol, (180/pi)*(xs_ltol(:, 3) + pi), 'Linewidth', 3)

figure('Color', 'w', 'Name', 'Top Linearized Cart Pendulum (Closed Loop)')
subplot(2, 1, 1), hold on, grid on, xlabel('Time [s]'), ylabel('Cart Position [m]'), title('Cart Position vs Time'), plot(ts_ltcl, xs_ltcl(:, 1), 'Linewidth', 3)
subplot(2, 1, 2), hold on, grid on, xlabel('Time [s]'), ylabel('Pendulum Angle [deg]'), title('Pendulum Angle vs Time'), plot(ts_ltcl, (180/pi)*(xs_ltcl(:, 3) + pi), 'Linewidth', 3)


