function [ phi, phi_t ] = DeterminingROA(ProjectOptions, phi0)
% For autonomous systems, because there are no continuous inputs, the reachable set calculation can be performed with convection instead of general HJ PDEs.

% Parameters:
%   projectoptions  Project Options for ROAtoolbox
%   phi0        Implicit surface function for initial set.
%   phi        Implicit surface function for ROA.
%   phi_t      phi for each period (include phi0).

% Retrieve the accuracy setting.
accuracy = ProjectOptions.accuracy;


%% Step 1: Generate the Grid.

% Create the grid.
g = ProjectOptions.Grid;


%% Step 2: Generate the Phi0 Initial Condition.

% Determine how to handle default input arguments.
if nargin < 2                                   % If the number of input arguments is less than two...
    
    % Retrieve the initial condition.
    phi0 = ProjectOptions.InitialCondition;
    
end


%% Step 3: Generate the Flow Field.

% Retrieve the flow field.
VecField = ProjectOptions.VectorField ;

% Reverse the flow field.
for ii = 1:size(VecField)
    VecField{ii} = -VecField{ii};
end


%% Step 4: Compute the Region of Attraction.

% Retrieve integration parameters.
tMax = ProjectOptions.IntegratorTime;               % End time for integrate.
tPlot = ProjectOptions.IntegratorShowTime;          % Period at which intermediate plots should be produced.

% A figure to put intermediate plots.
fig = figure;

% Visualize the  ROA.
switch(g.dim)
    case 1
        
        % Set the display type.
        displayType = 'plot';
        
    case 2
        
        % Set the display type.
        displayType = 'contour';
        
    case 3
        
        % Set the display type.
        displayType = 'surface';
        
    case 4
        
        % Set the display type.
        displayType = 'contour';
        
    otherwise
        
        % Throw an error.
        error('Can not display system with dimention: %s!', g.dim);
        
end

% Set the growth option (should generally be 1).
growOnly = 1;

% Compute the ROA.
[ phi, executionTime, phi_t ] = findReachSet(g, phi0, VecField, accuracy, tMax, tPlot, fig, growOnly, [], displayType);


%% ---------------------------------------------------------------------------
% Step 5, optional
% if g.dim == 2
%     figure;
%     level = [ 0 0 ];
%     contourf(g.xs{1}, g.xs{2}, -reach, level);
%     axis equal;
%     axis(g.axis);
%     xlabel('x_1'); ylabel('x_2');
%     % Clip back the axis bounds slightly so that the distortion caused
%     %   by the artificial boundary doesn't show.
%     % clip = [ +5, -5, +5, -5 ];
%     % axis(g.axis + clip);
% end

%%
fprintf('Total execution time %g seconds \n', executionTime);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------------
% The function findReachSet is written by Ian Mitchell, modified by YUAN
function [ final, executionTime, phi_t ] = findReachSet(g, initial, velocity, accuracy, tMax, tPlot, fig, growOnly, avoid, displayType)
% findReachSet: Compute reach or reach-avoid sets.
% Parameters:
%
%   g             Grid structure on which data was computed.
%   initial       Array containing the implicit surface for the target set.
%   velocity      Cell vector containing the convective flow field.
%   accuracy      Controls the order of approximations.
%                   'low'         Use odeCFL1 and upwindFirstFirst.
%                   'medium'      Use odeCFL2 and upwindFirstENO2.
%                   'high'        Use odeCFL3 and upwindFirstENO3.
%                   'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   tMax          Integration time (integrations start at t0 = 0).
%   tPlot         Time at which to produce intermediate figures.
%   fig           Handle to figure in which to produce intermediate figures.
%   growOnly      Boolean, specifies whether to restrict the temporal
%                   derivative of the implicit surface function so that the
%                   reachable set only grows.
%   avoid         Array containing the implicit surface for the escape set.
%                   Set to [] if there is no escape set.
%
%   final         Array containing the implicit surface function of the
%                   reach or reach-avoid set.
%   executionTime Time required to compute final (from cputime, in seconds).

% Ian Mitchell, 4/18/04

%---------------------------------------------------------------------------
% What level set should we view?
% level = 0;
level = 1;

% Visualize the  2D reachable set.
%displayType = 'contour';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

t0 = 0;                              % Start time.

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.grid = g;
schemeData.velocity = velocity;

%---------------------------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
    case 'low'
        schemeData.derivFunc = @upwindFirstFirst;
        integratorFunc = @odeCFL1;
    case 'medium'
        schemeData.derivFunc = @upwindFirstENO2;
        integratorFunc = @odeCFL2;
    case 'high'
        schemeData.derivFunc = @upwindFirstENO3;
        integratorFunc = @odeCFL3;
    case 'veryHigh'
        schemeData.derivFunc = @upwindFirstWENO5;
        integratorFunc = @odeCFL3;
    otherwise
        error('Unknown accuracy level %s', accuracy);
end

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.

if(growOnly)
    innerFunc = schemeFunc;
    innerData = schemeData;
    clear schemeFunc schemeData;
    
    % Wrap the true Hamiltonian inside the term approximation restriction.
    schemeFunc = @termRestrictUpdate;
    schemeData.innerFunc = innerFunc;
    schemeData.innerData = innerData;
    schemeData.positive = 0;
end

%---------------------------------------------------------------------------
if(isempty(avoid))
    data = initial;
else
    % Ensure that the initial data satisfies the avoid set.
    data = max(initial, avoid);
    
    % Set up data required for masking by the avoid set.
    %   Mask will be compared to vector form of data array used by integrator.
    schemeData.maskData = avoid(:);
    schemeData.maskFunc = @max;
    
    % Let the integrator know what function to call.
    integratorOptions = odeCFLset(integratorOptions, 'postTimestep', @postTimestepMask);
    
end

%---------------------------------------------------------------------------
% Initialize Display
tstr = [ 't = ' num2str(t0) ];
figure(fig);
h = visualizeLevelSet(g, data, displayType, level, tstr);
hold on;
%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
phi_t = cell(floor(tMax / tPlot) +1 , 1);
phi_t{1} = data;
iind = 2;
while(tMax - tNow > small * tMax)
    
    % Reshape data array into column vector for ode solver call.
    y0 = data(:);
    
    % How far to step?
    tSpan = [ tNow, min(tMax, tNow + tPlot) ];
    
    % Take a timestep.
    [ t, y ] = feval(integratorFunc, schemeFunc, tSpan, y0, integratorOptions, schemeData);
    tNow = t(end);
    
    % Get back the correctly shaped data array
    data = reshape(y, g.shape);
    
    % store phi for this period
    phi_t{iind} = data;
    iind = iind +1;
    save('phi_t', 'phi_t');
    
    if(pauseAfterPlot)
        % Wait for last plot to be digested.
        pause;
    end
    
    % Get correct figure.
    figure(fig);
    
    % Delete last visualization if necessary.
    if(deleteLastPlot)
        delete(h);
    end
    
    % Create new visualization.
    h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
    
end


executionTime = cputime - startTime;
fprintf('Execution time %g seconds\n', executionTime);

final = data;
