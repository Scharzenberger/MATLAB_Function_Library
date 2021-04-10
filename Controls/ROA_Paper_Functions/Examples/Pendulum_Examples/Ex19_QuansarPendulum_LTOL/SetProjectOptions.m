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
    
    % Retrieve the system inputs.
    [x1, x2, x3, x4] = xs{:};
    
    % Define the pendulum parameters.
    m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L1 = 1; L2 = 1; g = 9.81; etag = 1; etam = 1; kt = 1; km = 1; Kg = 1; Rm = 1; Vm = 0;
    
    % Compute the system derivatives (zero point is at the top).
    x1dot = x2;
    x2dot = ( -(12*k1)/(L1^2*m1 + 3*L1^2*m2) )*x1 + ( -(3*(8*L2*etag*etam*km*kt*Kg^2 + 8*L2*Rm*c1))/(2*L2*Rm*(L1^2*m1 + 3*L1^2*m2)) )*x2 + ( (3*(12*L1*Rm*k2 - 6*L1*L2*Rm*g*m2))/(2*L2*Rm*(L1^2*m1 + 3*L1^2*m2)) )*x3 + ( (18*L1*c2)/(L2*(L1^2*m1 + 3*L1^2*m2)) )*x4;
    x3dot = x4;
    x4dot = ( (18*L1*k1)/(L2*(L1^2*m1 + 3*L1^2*m2)) )*x1 + ( (3*(24*L1*L2*etag*etam*km*kt*m2*Kg^2 + 24*L1*L2*Rm*c1*m2))/(4*L2^2*Rm*m2*(L1^2*m1 + 3*L1^2*m2)) )*x2 + ( -(3*(4*L1^2*Rm*k2*m1 + 48*L1^2*Rm*k2*m2 - 24*L1^2*L2*Rm*g*m2^2 - 2*L1^2*L2*Rm*g*m1*m2))/(4*L2^2*Rm*m2*(L1^2*m1 + 3*L1^2*m2)) )*x3 + ( -(3*(4*L1^2*Rm*c2*m1 + 48*L1^2*Rm*c2*m2))/(4*L2^2*Rm*m2*(L1^2*m1 + 3*L1^2*m2)) )*x4;

    % Store the system derivatives.
    VectorField = { x1dot; x2dot; x3dot; x4dot };
    
    % -----------------------------------------------------------------------    


else   % scalar
    % Do not change
    % recursion
    VectorField = cell2mat (GenVecField( num2cell(xs) ));
end