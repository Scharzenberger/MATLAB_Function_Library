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
    
%     % -----------------------------------------------------------------------
% 
%     % Retrieve the system inputs.
%     [x1, x2, x3, x4, x5] = xs{:};
% 
%     % Define the pendulum parameters.
%     m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L1 = 1; L2 = 1; g = 9.81; etag = 1; etam = 1; kt = 1; km = 1; Kg = 1; Rm = 1; Vm = 0; xd = 0;
%     
%     % Define the state space controller.
%     K = [3.5302, 5.4000, -23.2219, -4.9432, -70.3667];      % Bottom Controller (tsettle = 10s).
% %     K = [-1252918.5164754043798894, 195.7499990213444221, -1354861.1859419560059905, -140876.7088980348780751, -7036672.3765125349164009];          % Bottom Controller (tsettle = 1s).
%     
%     % Compute the system derivatives.
%     f1 = x2;
%     f2 = ( -(12*k1)/(L1^2*m1 + 3*L1^2*m2) )*x1 + ( -(3*(8*L2*etag*etam*km*kt*Kg^2 + 8*L2*Rm*c1))/(2*L2*Rm*(L1^2*m1 + 3*L1^2*m2)) )*x2 + ( -(3*(12*L1*Rm*k2 + 6*L1*L2*Rm*g*m2))/(2*L2*Rm*(L1^2*m1 + 3*L1^2*m2)) )*x3 + ( -(18*L1*c2)/(L2*(L1^2*m1 + 3*L1^2*m2)) )*x4;
%     f3 = x4;
%     f4 = ( -(18*L1*k1)/(L2*(L1^2*m1 + 3*L1^2*m2)) )*x1 + ( -(3*(24*L1*L2*etag*etam*km*kt*m2*Kg^2 + 24*L1*L2*Rm*c1*m2))/(4*L2^2*Rm*m2*(L1^2*m1 + 3*L1^2*m2)) )*x2 + ( -(3*(4*L1^2*Rm*k2*m1 + 48*L1^2*Rm*k2*m2 + 24*L1^2*L2*Rm*g*m2^2 + 2*L1^2*L2*Rm*g*m1*m2))/(4*L2^2*Rm*m2*(L1^2*m1 + 3*L1^2*m2)) )*x3 + ( -(3*(4*L1^2*Rm*c2*m1 + 48*L1^2*Rm*c2*m2))/(4*L2^2*Rm*m2*(L1^2*m1 + 3*L1^2*m2)) )*x4;
% 
%     % Compute the closed loop system derivatives.
%     x1dot = f1;
%     x2dot = f2 + K(5)*x5 - K(1)*x1 - K(2)*x2 - K(3)*x3 - K(4)*x4;
%     x3dot = f3;
%     x4dot = f4;
%     x5dot = xd - x3;
% 
%     % Store the system derivatives.
%     VectorField = { x1dot; x2dot; x3dot; x4dot; x5dot };
%     
%     % -----------------------------------------------------------------------
    
        % Retrieve the system inputs.
    [x1, x2, x3, x4, x5] = xs{:};
    
    xd = 0;
    
    % Compute the closed loop system derivatives.
    x1dot = x2;
    x2dot = -6.5302*x1 + -11.4000*x2 + 1.1494*x3 + 0.4432*x4 + -70.3667*x5;
    x3dot = x4;
    x4dot = -4.5000*x1 + -9.0000*x2 + -47.8237*x3 + -9.7500*x4;
    x5dot = xd - x3;

    % Store the system derivatives.
    VectorField = { x1dot; x2dot; x3dot; x4dot; x5dot };
    
    
else   % scalar

    VectorField = cell2mat (GenVecField( num2cell(xs) ));
    
end
