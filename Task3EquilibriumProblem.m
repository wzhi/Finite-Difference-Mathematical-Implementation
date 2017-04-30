%% Solve Diffusion problem with Neumann- & Dirichlet BCs

%%
clear
FDLabFolders

%%
ToDo('Please work from DemoScripts/FDIncrementalDemo - not here!')
ToDo('Implement and display FDSystemMatrix')
dbtype FDSystemMatrix

%%
ToDo('Implement and display FDSystemVector')
dbtype FDSystemVector

%% Problem data
rho = 2; % reaction coefficient
kappa = [1 3; 2 4]; % diffusion coefficient
[exactSolution, domainSource, NeumannSource] = ...
    DataSteadyPoly(kappa, rho); % source terms
%%
% A non-uniformly spaced grid, which differs in each direction
x = lglspace(0, 2, 15); % grid point x-coordinates
y = lglspace(0, 3, 20); % grid point y-coordinates
xyGrid = FDGrid(x, y);
%%
% Apply Dirichlet conditions at an *arbitrary* selection of grid points.
Dirichlet = FDDirichlet( ...
    xyGrid, exactSolution, ...
    xyGrid.Indices([ 1, round(end/2),  end]));

%%
% Solve discretized problem:
U = FDSolve(xyGrid, kappa, rho, domainSource, NeumannSource, Dirichlet);
%%
fprintf('|relative errors|: %e\n', ...
    norm(Compare(FDEvaluate(exactSolution, xyGrid), U), inf))
