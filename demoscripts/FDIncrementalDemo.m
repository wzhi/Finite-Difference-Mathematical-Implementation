%% Script for testing FDM routines

%% Suggested Steps
% 
% # Please *publish* this script and study its output: Use |doc| to
% investigate any unfamiliar functions and speak to a tutor if you are
% unsure of anything.
% # Gradually increase the complexity of the configuration so as to build
% the capabilities of |FDSystemMatrix| and |FDSystemVector|. A list of such
% extensions follows.
%

%% Extensions to test
%
% # Remove Dirichlet conditions on interior grid points
% # Remove Dirichlet conditions on successive boundaries
% # Different spacing in the x-direction and y-direction
% # More than 3 points in each direction
% # Non-uniform spacing in each direction
% # Matrix-valued diffusion coefficient
% 
% _Start Simple and Proceed Mindfully to maximise productivity_!

%%
clear
FDLabFolders 
tol = 0;

%% Problem data
% Note that we are using a scalar value of the difffusion coefficient.
% Later on, we'll provide support for 2x2 arrays.
rho = 2; % reaction coefficient
kappa = [1 3;2 4]; % diffusion coefficient
[exactSolution, domainSource, NeumannSource] = ...
    DataSteadyPolyScalarKappa(kappa, rho);

%% Finite difference grid 
% You may like to have a quick look at the implementation of |FDGrid|.
x = lglspace(0, 2, 15); % grid point x-coordinates
y = lglspace(0, 3, 20); % grid point y-coordinates
xyGrid = FDGrid(x, y) %#ok<NOPTS> % Encodes a regular 2D grid
%%
disp(xyGrid.X) % x-coordinates of grid points
%%
disp(xyGrid.Y) % y-coordinates of grid points
%%
disp(xyGrid.Indices) % indices/labels of grid points
%% Dirichlet conditions
% Remove points on each boundary in sequence to test your 
% implementation of the finite difference equations. 
Dirichlet = FDDirichlet( ...
    xyGrid, ... % Encodes grid layout
    exactSolution, ... % Specified values at Dirichlet points
    ... xyGrid.Indices(2 : end - 1, 2 : end - 1), ... % interior grid points
    ...xyGrid.Indices([1 end], :), ... % points on north & south boundaries
    xyGrid.Indices([ 1, round(end/2),  end]) ... % points on east & west boundaries
    )  %#ok<NOPTS> 

%% Assemble and solve discrete equations
[U, A, b] = FDSolve(xyGrid, kappa, rho, ...
    domainSource, NeumannSource, Dirichlet);

%% Discrepancy with prescribed solution
discrepancy = Compare(FDEvaluate(exactSolution, xyGrid), U);
fprintf('Max. |relative error|: %1.2e\n', norm(discrepancy, inf))

%% System matrices
% Compare the system matrix and right-hand side with the system you have
% worked out on paper if the error norm is large.
disp(full([A b]))

%% Spatial distribution of errors
% Perhaps of use/interest if the discrepancy is large.
disp(discrepancy)

%% Visualization
clf, set(gca, 'FontSize', FontSize)
mesh(xyGrid.X, xyGrid.Y, U, 'EdgeAlpha', 0.1), hold on
plot3(xyGrid.X, xyGrid.Y, U, 'o', ...
    'MarkerSize', 3, 'MarkerFaceColor', 'r')
set(gca, 'XTick', x, 'XTickLabel', [], 'YTick', y, 'YTickLabel', [])
xlabel('x'), ylabel('y'), zlabel('u(x, y)')
plot3(xyGrid.X, xyGrid.Y, FDEvaluate(exactSolution,xyGrid), 'o', ...
    'MarkerSize', 5)
