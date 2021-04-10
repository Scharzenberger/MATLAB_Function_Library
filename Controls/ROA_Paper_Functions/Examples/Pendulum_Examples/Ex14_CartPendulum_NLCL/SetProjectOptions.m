function ROAToolboxOptions = SetProjectOptions()


%% Create ROA Computational Grid.

% The dimension of the grid should equal to to number of  state variables
ROAToolboxOptions.GridDimension = 5;

% Define the ROA grid domain.
% ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi];
% ROAToolboxOptions.GridRange = 2*[-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi];

ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -2*pi, 2*pi; -2*pi, 2*pi; -4*pi, 4*pi; -2*pi, 2*pi];
% ROAToolboxOptions.GridRange = [-2*pi, 2*pi; -8*pi, 8*pi; -2*pi, 2*pi; -4*pi, 4*pi; -2*pi, 2*pi];

% Define the ROA grid step size.
% ROAToolboxOptions.GridCellSize = [0.2; 0.2; 0.2; 0.2];
% ROAToolboxOptions.GridCellSize = [pi; pi; 0.2; 0.2];

ROAToolboxOptions.GridCellSize = [pi; pi; 0.2; 0.2; 0.2];
% ROAToolboxOptions.GridCellSize = [1; 1; 0.2; 0.2; 0.2];


%% Create ROA Initial Condition.

% Set the initial condition hyperspheroid radius.
ROAToolboxOptions.InitCircleRad = 1;

% Set the initial condition hyperspheroid centeroid.
% ROAToolboxOptions.InitCircleCen = [ 0; 0; 0; 0 ];
ROAToolboxOptions.InitCircleCen = [ 0; 0; 0; 0; 0 ];


%% Parameters for integrating

% Define the total integration time.
ROAToolboxOptions.IntegratorTime = 5;

% Define how often ROA results should be shown.
ROAToolboxOptions.IntegratorShowTime = 1;


%% Set the Computational Accuracy.

% Set the computational accuracy.
%    one of:  'low', 'medium', 'high', 'veryHigh'
ROAToolboxOptions.accuracy = 'medium';


%% Define the System of Equations.

% Define the vector field operator.
ROAToolboxOptions.VectorFieldOperator = @(x) GenVecField( x );


%% Assign the ROA Toolbox Options.

% ----------do not change -----------------

ROAToolboxOptions.Grid = GenerateGrid(ROAToolboxOptions);
fx = ROAToolboxOptions.VectorFieldOperator;
ROAToolboxOptions.VectorField = fx(ROAToolboxOptions.Grid.xs);
initialRadius = ROAToolboxOptions.InitCircleRad;
initialCenter = ROAToolboxOptions.InitCircleCen;
phi0 = shapeSphere(ROAToolboxOptions.Grid, initialCenter, initialRadius);
ROAToolboxOptions.InitialCondition = phi0;

% --------------------------------------------


%% Define the Dynamical System of Interest.

function VectorField = GenVecField( xs )

% Determine how to process the vector field input.
if iscell(xs)                       % If the vector field input is a cell...
    
    % -----------------------------------------------------------------------
    
    % Retrieve the system inputs.
    [x1, x2, x3, x4, x5] = xs{:};
    
    % Make the zero point of the closed loop nonlinear cart pendulum face upward.
    x3 = x3 + pi;
    
    % Define the pendulum parameters.
    m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L = 1; F = 0; g = 9.81; kmotor = 1; xd = 0;
    %     m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L = 1; F = 0; g = 9.81; kmotor = 1; xd = pi/4;
    
    % Define the state space controller.
%     K = [3.6774097132690682, 18.3318181818451649, -1244.5153035800524322, -247.1203927990935370, -580.5254812531755988];        % Top Controller (tsettle = 10s).
    K = [34.3763586229696259, 18.3318181818259802, 478.2803688165961375, 151.5989864640213227, 580.5254812494409862];           % Bottom Controller (tsettle = 10s).
    
%     K = [1915596.3311236712615937, 208.6818284108001080, 3037746.6962169660255313, 3481671.2530231652781367, -58052559.8825215548276901];        % Top Controller (tsettle = 1s).
%     K = [-2390735.6156045463867486, 208.6818257983125591, 6284470.0735367881134152, 4413231.4083734853193164, 58052554.8896254897117615];        % Bottom Controller (tsettle = 1s).
    
    % Compute the open loop system derivatives.
    f1 = x2;
    f2 = (-(7*k1)./(7*m1 - 3*m2*cos(x3).^2 + 7*m2)).*x1 + (-(7*c1)./(7*m1 - 3*m2*cos(x3).^2 + 7*m2)).*x2 + (-(6*k2*cos(x3))./(L*(7*m1 - 3*m2*cos(x3).^2 + 7*m2))).*x3 + (-(7*m2*sin(x3)*L^2 + 12*c2*cos(x3))./(2*L*(7*m1 - 3*m2*cos(x3).^2 + 7*m2))).*x4 + (14*F*L + 6*L*g*m2*cos(x3).*sin(x3))./(2*L*(7*m1 - 3*m2*cos(x3).^2 + 7*m2));
    f3 = x4;
    f4 = (-(6*k1*cos(x3))./(L*(7*m1 - 3*m2*cos(x3).^2 + 7*m2))).*x1 + (-(6*c1*cos(x3))./(L*(7*m1 - 3*m2*cos(x3).^2 + 7*m2))).*x2 + (-(3*(4*k2*m1 + 4*k2*m2))./(L^2*m2*(7*m1 - 3*m2*cos(x3).^2 + 7*m2))).*x3 + (-(3*(cos(x3).*sin(x3)*L^2*m2^2 + 4*c2*m2 + 4*c2*m1))./(L^2*m2*(7*m1 - 3*m2*cos(x3).^2 + 7*m2))).*x4 + (3*(2*F*L*m2*cos(x3) + 2*L*g*m2^2*sin(x3) + 2*L*g*m1*m2*sin(x3)))./(L^2*m2*(7*m1 - 3*m2*cos(x3).^2 + 7*m2));
    
    % Compute the closed loop system derivatives.
    x1dot = f1;
    x2dot = f2 + kmotor*(K(5)*x5 - K(1)*x1 - K(2)*x2 - K(3)*x3 - K(4)*x4);
    x3dot = f3;
    x4dot = f4;
    x5dot = xd - x3;
    
    % Store the system derivatives.
    VectorField = { x1dot; x2dot; x3dot; x4dot; x5dot };

%     % Define the pendulum parameters.
%     m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L1 = 1; L2 = 1; g = 9.81; etag = 1; etam = 1; kt = 1; km = 1; Kg = 1; Rm = 1; Vm = 0; xd = 0;
%     
%     % Define the state space controller.
%     K = [1.5603, 5.4000, 96.9413, 15.7895, 70.3667];        % Top Controller (tsettle = 10s).
% 
%     % Compute the system derivatives.
%     f1 = x2;
%     f2 = -(3*(8*L2*Rm*c1*x2 + 8*L2*Rm*k1*x1 + 12*L1*Rm*c2*x4.*cos(x3) + 12*L1*Rm*k2*x3.*cos(x3) + 2*L2^3*Rm*m2*x2.*x4.*sin(2*x3) + 3*L1*L2*Rm*g*m2*sin(2*x3) - 8*Kg*L2*Vm*etag*etam*kt + 4*L1*L2^2*Rm*m2*x4.^2.*sin(x3) - 3*L1*L2^2*Rm*m2*x2.^2.*(sin(x3) - sin(x3).^3) + 8*Kg^2*L2*etag*etam*km*kt*x2))./(2*L2*Rm*(L1^2*m1 + 12*L1^2*m2 + 3*L2^2*m2 - 9*L1^2*m2*cos(x3).^2 - 3*L2^2*m2*cos(x3).^2));
%     f3 = x4;
%     f4 = -(3*(4*L1^2*Rm*c2*m1*x4 - (3*L2^4*Rm*m2^2*x2.^2.*sin(2*x3))/2 + 48*L1^2*Rm*c2*m2*x4 + 12*L2^2*Rm*c2*m2*x4 + 4*L1^2*Rm*k2*m1*x3 + 48*L1^2*Rm*k2*m2*x3 + 12*L2^2*Rm*k2*m2*x3 + 6*L2^3*Rm*g*m2^2*sin(x3).^3 + 3*L2^4*Rm*m2^2*x2.^2.*cos(x3).^3.*sin(x3) + 24*L1^2*L2*Rm*g*m2^2*sin(x3) - 12*L2^2*Rm*c2*m2*x4.*cos(x3).^2 - 12*L2^2*Rm*k2*m2*x3.*cos(x3).^2 - 6*L1^2*L2^2*Rm*m2^2*x2.^2.*sin(2*x3) + 6*L1^2*L2^2*Rm*m2^2*x4.^2.*sin(2*x3) - (L1^2*L2^2*Rm*m1*m2*x2.^2.*sin(2*x3))/2 + 12*L1*L2^3*Rm*m2^2*x2.*x4.*sin(x3) + 24*L1*L2*Rm*c1*m2*x2.*cos(x3) + 24*L1*L2*Rm*k1*m2*x1.*cos(x3) - 12*L1*L2^3*Rm*m2^2*x2.*x4.*sin(x3).^3 + 2*L1^2*L2*Rm*g*m1*m2*sin(x3) - 24*Kg*L1*L2*Vm*etag*etam*kt*m2*cos(x3) + 24*Kg^2*L1*L2*etag*etam*km*kt*m2*x2.*cos(x3)))./(4*L2^2*Rm*m2*(L1^2*m1 + 12*L1^2*m2 + 3*L2^2*m2 - 9*L1^2*m2*cos(x3).^2 - 3*L2^2*m2*cos(x3).^2));
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
    
    % -----------------------------------------------------------------------
    
    
    %     % -----------------------------------------------------------------------
    %
    %     % Retrieve the system inputs.
    %     [x1, x2, x3, x4] = xs{:};
    %
    %     % Make the zero point of the open loop nonlinear chart pendulum face downward.
    %     x3 = x3 + pi;
    %
    %     % Define the pendulum parameters.
    %     m1 = 1; m2 = 1; c1 = 1; c2 = 1; k1 = 1; k2 = 0; L = 1; F = 0; g = 9.81; kmotor = 1;
    %
    %     % Define the state space gain.
    %     K = [41730.7994636354560498, 150.1818181793833560, 173583.1544857693370432, -62465.6046231265427195];
    % %     K = [0, 0, 0, 0];
    %
    %     % Compute the system derivatives.
    %     x1dot = x2;
    %     x2dot = - (7*L.*m2.*sin(x3).*x4.^2)./(- 6*m2.*cos(x3).^2 + 14*m1 + 14*m2) - (6*c2.*cos(x3).*x4)./(L.*(- 3*m2.*cos(x3).^2 + 7*m1 + 7*m2)) + (14*F.*L + 6*L.*g.*m2.*cos(x3).*sin(x3))./(2*L.*(- 3*m2.*cos(x3).^2 + 7*m1 + 7*m2)) - (7*c1.*x2)./(- 3*m2.*cos(x3).^2 + 7*m1 + 7*m2) - (7*k1.*x1)./(- 3*m2.*cos(x3).^2 + 7*m1 + 7*m2) - (6*k2.*x3.*cos(x3))./(L.*(- 3*m2.*cos(x3).^2 + 7*m1 + 7*m2)) + kmotor*K(1)*x1 + kmotor*K(2)*x2 + kmotor*K(3)*x3 + kmotor*K(4)*x4;
    %     x3dot = x4;
    %     x4dot = (-(6*k1.*cos(x3))./(L.*(7*m1 - 3*m2.*cos(x3).^2 + 7*m2))).*x1 + (-(6*c1.*cos(x3))./(L.*(7*m1 - 3*m2.*cos(x3).^2 + 7*m2))).*x2 + (-(3*(4*k2.*m1 + 4*k2.*m2))./(L.^2.*m2.*(7*m1 - 3*m2.*cos(x3).^2 + 7*m2))).*x3 + (-(3*m2.*cos(x3).*sin(x3))./(7*m1 - 3*m2.*cos(x3).^2 + 7*m2)).*x4.^2 + (-(3*(4*c2.*m1 + 4*c2.*m2))./(L.^2.*m2.*(7*m1 - 3*m2.*cos(x3).^2 + 7*m2))).*x4 + (3*(2*F.*L.*m2.*cos(x3) + 2*L.*g.*m2.^2.*sin(x3) + 2*L.*g.*m1.*m2.*sin(x3)))./(L.^2.*m2.*(7*m1 - 3*m2.*cos(x3).^2 + 7*m2));
    %
    %     % Store the system derivatives.
    %     VectorField = { x1dot; x2dot; x3dot; x4dot };
    %
    %     % -----------------------------------------------------------------------
    
    
else
    
    VectorField = cell2mat( GenVecField( num2cell(xs) ) );
end
