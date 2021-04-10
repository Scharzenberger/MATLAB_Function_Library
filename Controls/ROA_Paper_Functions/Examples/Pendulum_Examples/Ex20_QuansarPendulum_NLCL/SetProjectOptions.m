function ROAToolboxOptions = SetProjectOptions()


%% Parameters for creating computational grid

% The dimension of the grid should equal to to number of  state variables
ROAToolboxOptions.GridDimension = 5;

% Define the grid domain.
ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -4*pi, 4*pi; -2*pi, 2*pi; -4*pi, 4*pi; -2*pi, 2*pi];

% Define the grid step size.
% ROAToolboxOptions.GridCellSize = [0.2; 0.2; 0.2; 0.2];
ROAToolboxOptions.GridCellSize = [pi; pi; 0.2; 0.2; 0.2];


%% Initial condition

% Define the radius of the initial condition hyperspheriod.
ROAToolboxOptions.InitCircleRad = 1;

% Define the center of the initial condition hyperspheriod.
ROAToolboxOptions.InitCircleCen = [ 0; 0; 0; 0; 0 ];


%% Parameters for integrating

% Define the intergration time.
ROAToolboxOptions.IntegratorTime = 5;

% Define the interval at which to plot intermediate results.
ROAToolboxOptions.IntegratorShowTime = 1;


%% Set computational accuracy

% Define the accuracy of the numerical method.
ROAToolboxOptions.accuracy = 'medium';


%% System equation

% Define the system vector field.
ROAToolboxOptions.VectorFieldOperator = @(x) GenVecField( x );


%% For efficient, do not change
% ----------do not change -----------------
ROAToolboxOptions.Grid = GenerateGrid(ROAToolboxOptions);
fx = ROAToolboxOptions.VectorFieldOperator;
ROAToolboxOptions.VectorField = fx(ROAToolboxOptions.Grid.xs);
initialRadius = ROAToolboxOptions.InitCircleRad;
initialCenter = ROAToolboxOptions.InitCircleCen;
phi0 = shapeSphere(ROAToolboxOptions.Grid, initialCenter, initialRadius);
ROAToolboxOptions.InitialCondition = phi0;
% --------------------------------------------


%% System equation,
function VectorField = GenVecField( xs )
% ---------Modify according to the undetlying system------------
% Compute the vector field of the system.
% This function should be element-wise.
% Parameters:
% xs is either a cell(n*1) contain n members each of which is a
%   multidimansional array, or a n*1 array containing the state variable.
% VectorField is the same type as xs.
%
% YUAN Guoqiang, Oct, 2016
%
if iscell(xs)  % element-wise
    
    % -----------------------------------------------------------------------
    
%     % Retrieve the system inputs.
%     [x1, x2, x3, x4] = xs{:};
%     
%     % Define the pendulum parameters.
%     m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L1 = 1; L2 = 1; g = 9.81; etag = 1; etam = 1; kt = 1; km = 1; Kg = 1; Rm = 1; Vm = 0;
%     
%     % Compute the system derivatives.
%     x1dot = x2;
%     x2dot = -(3*(8*L2*Rm*c1*x2 + 8*L2*Rm*k1*x1 + 12*L1*Rm*c2*x4.*cos(x3) + 12*L1*Rm*k2*x3.*cos(x3) + 2*L2^3*Rm*m2*x2.*x4.*sin(2*x3) + 3*L1*L2*Rm*g*m2*sin(2*x3) - 8*Kg*L2*Vm*etag*etam*kt + 4*L1*L2^2*Rm*m2*x4.^2.*sin(x3) - 3*L1*L2^2*Rm*m2*x2.^2.*(sin(x3) - sin(x3).^3) + 8*Kg^2*L2*etag*etam*km*kt*x2))./(2*L2*Rm*(L1^2*m1 + 12*L1^2*m2 + 3*L2^2*m2 - 9*L1^2*m2*cos(x3).^2 - 3*L2^2*m2*cos(x3).^2));
%     x3dot = x4;
%     x4dot = -(3*(4*L1^2*Rm*c2*m1*x4 - (3*L2^4*Rm*m2^2*x2.^2.*sin(2*x3))/2 + 48*L1^2*Rm*c2*m2*x4 + 12*L2^2*Rm*c2*m2*x4 + 4*L1^2*Rm*k2*m1*x3 + 48*L1^2*Rm*k2*m2*x3 + 12*L2^2*Rm*k2*m2*x3 + 6*L2^3*Rm*g*m2^2*sin(x3).^3 + 3*L2^4*Rm*m2^2*x2.^2.*cos(x3).^3.*sin(x3) + 24*L1^2*L2*Rm*g*m2^2*sin(x3) - 12*L2^2*Rm*c2*m2*x4.*cos(x3).^2 - 12*L2^2*Rm*k2*m2*x3.*cos(x3).^2 - 6*L1^2*L2^2*Rm*m2^2*x2.^2.*sin(2*x3) + 6*L1^2*L2^2*Rm*m2^2*x4.^2.*sin(2*x3) - (L1^2*L2^2*Rm*m1*m2*x2.^2.*sin(2*x3))/2 + 12*L1*L2^3*Rm*m2^2*x2.*x4.*sin(x3) + 24*L1*L2*Rm*c1*m2*x2.*cos(x3) + 24*L1*L2*Rm*k1*m2*x1.*cos(x3) - 12*L1*L2^3*Rm*m2^2*x2.*x4.*sin(x3).^3 + 2*L1^2*L2*Rm*g*m1*m2*sin(x3) - 24*Kg*L1*L2*Vm*etag*etam*kt*m2*cos(x3) + 24*Kg^2*L1*L2*etag*etam*km*kt*m2*x2.*cos(x3)))./(4*L2^2*Rm*m2*(L1^2*m1 + 12*L1^2*m2 + 3*L2^2*m2 - 9*L1^2*m2*cos(x3).^2 - 3*L2^2*m2*cos(x3).^2));
%     
%     % Store the system derivatives.
%     VectorField = { x1dot; x2dot; x3dot; x4dot};
    
    % Retrieve the system inputs.
    [x1, x2, x3, x4, x5] = xs{:};

    % Define the pendulum parameters.
    m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L1 = 1; L2 = 1; g = 9.81; etag = 1; etam = 1; kt = 1; km = 1; Kg = 1; Rm = 1; Vm = 0; xd = 0;
    
    % Define the state space controller.
    K = [3.5302, 5.4000, -23.2219, -4.9432, -70.3667];      % Bottom Controller (tsettle = 10s).
%     K = [1.5603, 5.4000, 96.9413, 15.7895, 70.3667];        % Top Controller (tsettle = 10s).
    
    % Compute the system derivatives.
    f1 = x2;
    f2 = -(3*(8*L2*Rm*c1*x2 + 8*L2*Rm*k1*x1 + 12*L1*Rm*c2*x4.*cos(x3) + 12*L1*Rm*k2*x3.*cos(x3) + 2*L2^3*Rm*m2*x2.*x4.*sin(2*x3) + 3*L1*L2*Rm*g*m2*sin(2*x3) - 8*Kg*L2*Vm*etag*etam*kt + 4*L1*L2^2*Rm*m2*x4.^2.*sin(x3) - 3*L1*L2^2*Rm*m2*x2.^2.*(sin(x3) - sin(x3).^3) + 8*Kg^2*L2*etag*etam*km*kt*x2))./(2*L2*Rm*(L1^2*m1 + 12*L1^2*m2 + 3*L2^2*m2 - 9*L1^2*m2*cos(x3).^2 - 3*L2^2*m2*cos(x3).^2));
    f3 = x4;
    f4 = -(3*(4*L1^2*Rm*c2*m1*x4 - (3*L2^4*Rm*m2^2*x2.^2.*sin(2*x3))/2 + 48*L1^2*Rm*c2*m2*x4 + 12*L2^2*Rm*c2*m2*x4 + 4*L1^2*Rm*k2*m1*x3 + 48*L1^2*Rm*k2*m2*x3 + 12*L2^2*Rm*k2*m2*x3 + 6*L2^3*Rm*g*m2^2*sin(x3).^3 + 3*L2^4*Rm*m2^2*x2.^2.*cos(x3).^3.*sin(x3) + 24*L1^2*L2*Rm*g*m2^2*sin(x3) - 12*L2^2*Rm*c2*m2*x4.*cos(x3).^2 - 12*L2^2*Rm*k2*m2*x3.*cos(x3).^2 - 6*L1^2*L2^2*Rm*m2^2*x2.^2.*sin(2*x3) + 6*L1^2*L2^2*Rm*m2^2*x4.^2.*sin(2*x3) - (L1^2*L2^2*Rm*m1*m2*x2.^2.*sin(2*x3))/2 + 12*L1*L2^3*Rm*m2^2*x2.*x4.*sin(x3) + 24*L1*L2*Rm*c1*m2*x2.*cos(x3) + 24*L1*L2*Rm*k1*m2*x1.*cos(x3) - 12*L1*L2^3*Rm*m2^2*x2.*x4.*sin(x3).^3 + 2*L1^2*L2*Rm*g*m1*m2*sin(x3) - 24*Kg*L1*L2*Vm*etag*etam*kt*m2*cos(x3) + 24*Kg^2*L1*L2*etag*etam*km*kt*m2*x2.*cos(x3)))./(4*L2^2*Rm*m2*(L1^2*m1 + 12*L1^2*m2 + 3*L2^2*m2 - 9*L1^2*m2*cos(x3).^2 - 3*L2^2*m2*cos(x3).^2));
    
    % Compute the closed loop system derivatives.
    x1dot = f1;
    x2dot = f2 + K(5)*x5 - K(1)*x1 - K(2)*x2 - K(3)*x3 - K(4)*x4;
    x3dot = f3;
    x4dot = f4;
    x5dot = xd - x3;

    % Store the system derivatives.
    VectorField = { x1dot; x2dot; x3dot; x4dot; x5dot };

    % -----------------------------------------------------------------------
else   % scalar
    % Do not change
    % recursion
    VectorField = cell2mat (GenVecField( num2cell(xs) ));
end
