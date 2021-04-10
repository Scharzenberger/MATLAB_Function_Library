function ROAToolboxOptions = SetProjectOptions()
% Setting the Project Options
% Edit the file to change the options for the underlying system
%
% YUAN Guoqiang, Oct, 2016
%
%% Parameters for creating computational grid
% The dimension of the grid should equal to to number of  state variables
ROAToolboxOptions.GridDimension = 2;
% The toolboxLS-1.1.1 work best if all the dimension in the problem are
%   approximately the same size: for example, the grid ranges and cell
%   widths should be within an order of magnitude of one another. So the
%   next two parameters should be approximately same.
% The range of each dimension i.e. [x1min, x1max; x2min, x2max; x3min, x3max ...]
ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -2*pi, 2*pi];

% A column vector specifying the size of grid cell in each dimension
%  i.e. [x1cellsize; x2cellsize; x3cellsize ...]
ROAToolboxOptions.GridCellSize = [0.2; 0.2];
% A column vector specifying the number of grid nodes in each dimension.
%   this parameter can be automatically generated from GridRange and
%   GridCellSize.
% ROAToolboxOptions.GridNodeNumber = [101; 101];

%% Initial condition
% In the version we use a circle as the initial condition because it always work well enough.
%   The initial circle is an initial estimate of the ROA
%   InitCircleRad is the  Radius of initial circle (positive).
%   InitCircleCen is the  center of initial circle (column vector [x1; x2; x3 ...]).
ROAToolboxOptions.InitCircleRad = 1;
ROAToolboxOptions.InitCircleCen = [ 0; 0 ];

%% Parameters for integrating
% End time for integrate.
ROAToolboxOptions.IntegratorTime = 5;
% Period at which intermediate results should be produced.
ROAToolboxOptions.IntegratorShowTime = 1;

%% Set computational accuracy
% accuracy     Controls the order of approximations.
%    one of:  'low', 'medium', 'high', 'veryHigh'
ROAToolboxOptions.accuracy = 'medium';

%% System equation
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
    
    % Retrieve the components of the input vector.
    [x1, x2] = xs{:};
    
    % Define the pendulum parameters.
    %     m = 1; c = 1; k = 1; L = 1; g = 9.81; kmotor = 1;
    %     m = 1; c = 1; k = 0; L = 1; g = 9.81; kmotor = 1; kp = 1;
    
    %     m = 1; c = 1; k = 0; L = 1; g = 9.81; kmotor = 1; kp = 20;
    m = 1; c = 1; k = 0; L = 1; g = 9.81; kmotor = 1; kp = 2.0;
    %     m = 1; c = 1; k = 0; L = 1; g = 9.81; kmotor = 1; kp = 0;
    
    
    % Set the desired pendulum angle (zero is at the top).
    thetad = 0;
    
    % Define the nonlinear simple pendulum dynamics (zero is at top).
    x1dot = x2;
    x2dot = ((m*g*L - (k + kmotor*kp))/(m*(L^2)))*x1 + (-c/(m*(L^2)))*x2 + (kmotor*kp/(m*(L^2)))*thetad;
    
    VectorField = { x1dot;   x2dot};
    
    
    % % State Space Controller
    
    %     % Retrieve the components of the input vector.
    %     [x1, x2, x3] = xs{:};
    %
    % %     % Define the pendulum parameters.
    % %     m = 1; L = 1; g = 9.81; c = 1; ktau = 0;
    %
    %     % Define the open loop simple pendulum linearized at bottom (Example 6) (Zero is at bottom).
    %     x1dot = x2;
    %     x2dot = -2652.75*x1 + -99*x2 + 10023.75*x3;
    %     x3dot = -x1;
    %
    %     % Define the vector field.
    %     VectorField = { x1dot; x2dot; x3dot };
    
    % -----------------------------------------------------------------------
else   % scalar
    % Do not change
    % recursion
    VectorField = cell2mat (GenVecField( num2cell(xs) ));
end
