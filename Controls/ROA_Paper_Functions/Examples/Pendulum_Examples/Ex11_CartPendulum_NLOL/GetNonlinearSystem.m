%% Cart Pendulum: Get Nonlinear System

% This script computes the nonlinear system of equations that describe the nonlinear cart pendulum system.

% Clear Everything.
clear, close('all'), clc


%% Define the Cart Pendulum System of Equations

% Define the symbolic variables.
syms m1 m2 c1 c2 k1 k2 L g F x1 x2 x3 x4 x1dot x2dot x3dot x4dot

% Define the system of equations.
eq1 = (m1 + m2)*x2dot + c1*x2 + (1/2)*m2*L*x4*sin(x3) - (1/2)*m2*L*x4dot*cos(x3) == F;
eq2 = ((7/12)*m2*(L^2))*x4dot + c2*x4 - (1/2)*m2*g*L*sin(x3) - (1/2)*m2*L*cos(x3)*x2dot == 0;

% eq1 = (m1 + m2)*x2dot + c1*x2 - (1/2)*m2*L*x4dot == F;
% eq2 = ((7/12)*m2*(L^2))*x4dot + c2*x4 - (1/2)*m2*g*L*x3 - (1/2)*m2*L*x2dot == 0;

% Solve the system of equations.
sol = solve([eq1, eq2], [x2dot, x4dot]);

% Retrieve the solved system of equations.
x1dot = x2;
x2dot = sol.x2dot;
x3dot = x4;
x4dot = sol.x4dot;


% %% Compute the System Jacobian.
% 
% % Set the equilibrium point.
% xeq_bottom = [0; 0; pi; 0];
% xeq_top = [0; 0; 0; 0];
% 
% % Compute the system jacobian.
% J = [diff(x1dot, x1) diff(x1dot, x2) diff(x1dot, x3) diff(x1dot, x4); diff(x2dot, x1) diff(x2dot, x2) diff(x2dot, x3) diff(x2dot, x4); diff(x3dot, x1) diff(x3dot, x2) diff(x3dot, x3) diff(x3dot, x4); diff(x4dot, x1) diff(x4dot, x2) diff(x4dot, x3) diff(x4dot, x4)];
% 
% % Evaluate the jacobian at the equilibrium point.
% J_bottom = subs(J, [x1, x2, x3, x4], xeq_bottom');
% 
% % Evaluate the jacobian at the equilibrium point.
% J_top = subs(J, [x1, x2, x3, x4], xeq_top');

