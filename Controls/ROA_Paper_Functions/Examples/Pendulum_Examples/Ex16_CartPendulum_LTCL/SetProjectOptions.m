function ROAToolboxOptions = SetProjectOptions()


%% Define the Computational Grid.

% The dimension of the grid should equal to to number of  state variables
ROAToolboxOptions.GridDimension = 5;

% Define the grid domains.
% ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi];
% ROAToolboxOptions.GridRange = 2*[-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi];

ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -4*pi, 4*pi; -2*pi, 2*pi];


% Define the grid step size.
% ROAToolboxOptions.GridCellSize = [pi; pi; 0.2; 0.2];
ROAToolboxOptions.GridCellSize = [pi; pi; 0.2; 0.2; 0.2];


%% Define the Initial Condition.

% Define the radus of the initial condition spheroid.
ROAToolboxOptions.InitCircleRad = 1;

% Define the center of the initial condition spheroid.
ROAToolboxOptions.InitCircleCen = [ 0; 0; 0; 0; 0 ];


%% Define Integration Parameters.

% End time for integrate.
ROAToolboxOptions.IntegratorTime = 5;

% Period at which intermediate results should be produced.
ROAToolboxOptions.IntegratorShowTime = 1;


%% Define Computational Accuracy.

% accuracy     Controls the order of approximations.
%    one of:  'low', 'medium', 'high', 'veryHigh'
ROAToolboxOptions.accuracy = 'medium';


%% Define the System Dynamics.

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

if iscell(xs)
    
    
%     % ---------- SS Controller ----------
%     
%     % Retrieve the system inputs.
%     [x1, x2, x3, x4] = xs{:};
%     
%     % Define the pendulum parameters.
%     m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 0; k2 = 0; L = 1; g = 9.81; kmotor = 1; xd = 0;
% 
%     % Define the state space gain.
%     K = [-50578.5550492562106228, 150.1818181868150930, -486269.7224879938876256, -106808.6186516829475295];
%     
%     % Define the open loop system matrices.
%     A = [0, 1, 0, 0; -(7*k1)/(7*m1 + 4*m2), -(7*c1)/(7*m1 + 4*m2), -(12*k2 - 6*L*g*m2)/(2*L*(7*m1 + 4*m2)), -(6*c2)/(L*(7*m1 + 4*m2)); 0, 0, 0, 1; -(6*k1)/(L*(7*m1 + 4*m2)), -(6*c1)/(L*(7*m1 + 4*m2)), -(3*(4*k2*m1 + 4*k2*m2 - 2*L*g*m2^2 - 2*L*g*m1*m2))/(L^2*m2*(7*m1 + 4*m2)), -(3*(4*c2*m1 + 4*c2*m2))/(L^2*m2*(7*m1 + 4*m2))];
%     B = [0; kmotor; 0; 0];
% 
%     % Define closed loop system matrices.
%     Acl = A - B*K;
%     Bcl = B;
%     
%     % Compute the system derivatives (zero point is at the bottom).
%     x1dot = Acl(1, 1)*x1 + Acl(1, 2)*x2 + Acl(1, 3)*x3 + Acl(1, 4)*x4 + Bcl(1)*xd;
%     x2dot = Acl(2, 1)*x1 + Acl(2, 2)*x2 + Acl(2, 3)*x3 + Acl(2, 4)*x4 + Bcl(2)*xd;
%     x3dot = Acl(3, 1)*x1 + Acl(3, 2)*x2 + Acl(3, 3)*x3 + Acl(3, 4)*x4 + Bcl(3)*xd;
%     x4dot = Acl(4, 1)*x1 + Acl(4, 2)*x2 + Acl(4, 3)*x3 + Acl(4, 4)*x4 + Bcl(4)*xd;
% 
%     % Store the system derivatives.
%     VectorField = { x1dot; x2dot; x3dot; x4dot };
    

    % ---------- SSI Controller ----------
    
    % Retrieve the system inputs.
    [x1, x2, x3, x4, x5] = xs{:};
    
    % Define the pendulum parameters.
    m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L = 1; g = 9.81; kmotor = 1; xd = 0;

    % Define the state space gain.
    K = [3.6774097132690682, 18.3318181818451649, -1244.5153035800524322, -247.1203927990935370, -580.5254812531755988];

    % Define the open loop system matrices.
    A = [0, 1, 0, 0; -(7*k1)/(7*m1 + 4*m2), -(7*c1)/(7*m1 + 4*m2), -(12*k2 - 6*L*g*m2)/(2*L*(7*m1 + 4*m2)), -(6*c2)/(L*(7*m1 + 4*m2)); 0, 0, 0, 1; -(6*k1)/(L*(7*m1 + 4*m2)), -(6*c1)/(L*(7*m1 + 4*m2)), -(3*(4*k2*m1 + 4*k2*m2 - 2*L*g*m2^2 - 2*L*g*m1*m2))/(L^2*m2*(7*m1 + 4*m2)), -(3*(4*c2*m1 + 4*c2*m2))/(L^2*m2*(7*m1 + 4*m2))];
    B = [0; kmotor; 0; 0];
    C = [0 0 1 0];
    
    % Define closed loop system matrices.
    Acl = [A - B*K(1:end-1), B*K(end); -C, 0];
    Bcl = [zeros(size(A, 1), 1); 1];
    
    % Compute the system derivatives (zero point is at the bottom).
    x1dot = Acl(1, 1)*x1 + Acl(1, 2)*x2 + Acl(1, 3)*x3 + Acl(1, 4)*x4 + Acl(1, 5)*x5 + Bcl(1)*xd;
    x2dot = Acl(2, 1)*x1 + Acl(2, 2)*x2 + Acl(2, 3)*x3 + Acl(2, 4)*x4 + Acl(2, 5)*x5 + Bcl(2)*xd;
    x3dot = Acl(3, 1)*x1 + Acl(3, 2)*x2 + Acl(3, 3)*x3 + Acl(3, 4)*x4 + Acl(3, 5)*x5 + Bcl(3)*xd;
    x4dot = Acl(4, 1)*x1 + Acl(4, 2)*x2 + Acl(4, 3)*x3 + Acl(4, 4)*x4 + Acl(4, 5)*x5 + Bcl(4)*xd;
    x5dot = Acl(5, 1)*x1 + Acl(5, 2)*x2 + Acl(5, 3)*x3 + Acl(5, 4)*x4 + Acl(5, 5)*x5 + Bcl(5)*xd;

    % Store the system derivatives.
    VectorField = { x1dot; x2dot; x3dot; x4dot; x5dot };
    
else

    
    VectorField = cell2mat (GenVecField( num2cell(xs) ));
    
end
