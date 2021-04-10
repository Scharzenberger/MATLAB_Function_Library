function ROAToolboxOptions = SetProjectOptions()
% Setting the Project Options
% Edit the file to change the options for the underlying system
%
% YUAN Guoqiang, Oct, 2016
%
%% Parameters for creating computational grid

% The dimension of the grid should equal to to number of  state variables
ROAToolboxOptions.GridDimension = 4;

% The toolboxLS-1.1.1 work best if all the dimension in the problem are
%   approximately the same size: for example, the grid ranges and cell
%   widths should be within an order of magnitude of one another. So the
%   next two parameters should be approximately same.
% The range of each dimension i.e. [x1min, x1max; x2min, x2max; x3min, x3max ...]
% ROAToolboxOptions.GridRange = [-10, 10; -10, 10; -10, 10; -10, 10];
ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi];
% ROAToolboxOptions.GridRange = 2*[-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi];

% A column vector specifying the size of grid cell in each dimension
%  i.e. [x1cellsize; x2cellsize; x3cellsize ...]
% ROAToolboxOptions.GridCellSize = [0.2; 0.2; 0.2; 0.2];
ROAToolboxOptions.GridCellSize = [pi; pi; 0.2; 0.2];


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
ROAToolboxOptions.InitCircleCen = [ 0; 0; 0; 0 ];

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

if iscell(xs)  % element-wise
    
    % -----------------------------------------------------------------------
    
    % Retrieve the system inputs.
    [x1, x2, x3, x4] = xs{:};
    
    % Define the pendulum parameters.
    m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L = 1; g = 9.81; kmotor = 1; xd = 0;
    
    % Define the open loop system matrices.
    A = [0, 1, 0, 0; -(7*k1)/(7*m1 + 4*m2), -(7*c1)/(7*m1 + 4*m2), -(12*k2 - 6*L*g*m2)/(2*L*(7*m1 + 4*m2)), -(6*c2)/(L*(7*m1 + 4*m2)); 0, 0, 0, 1; -(6*k1)/(L*(7*m1 + 4*m2)), -(6*c1)/(L*(7*m1 + 4*m2)), -(3*(4*k2*m1 + 4*k2*m2 - 2*L*g*m2^2 - 2*L*g*m1*m2))/(L^2*m2*(7*m1 + 4*m2)), -(3*(4*c2*m1 + 4*c2*m2))/(L^2*m2*(7*m1 + 4*m2))];
    B = [0; kmotor; 0; 0];
    
    % Compute the system derivatives (zero point is at the bottom).
    x1dot = A(1, 1)*x1 + A(1, 2)*x2 + A(1, 3)*x3 + A(1, 4)*x4 + B(1)*xd;
    x2dot = A(2, 1)*x1 + A(2, 2)*x2 + A(2, 3)*x3 + A(2, 4)*x4 + B(2)*xd;
    x3dot = A(3, 1)*x1 + A(3, 2)*x2 + A(3, 3)*x3 + A(3, 4)*x4 + B(3)*xd;
    x4dot = A(4, 1)*x1 + A(4, 2)*x2 + A(4, 3)*x3 + A(4, 4)*x4 + B(4)*xd;
    
    % Store the system derivatives.
    VectorField = { x1dot; x2dot; x3dot; x4dot };
    
    % -----------------------------------------------------------------------
    
    
else   % scalar

    VectorField = cell2mat (GenVecField( num2cell(xs) ));
    
end
