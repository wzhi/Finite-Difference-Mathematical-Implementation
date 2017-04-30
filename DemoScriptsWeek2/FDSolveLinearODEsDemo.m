%% The theta method

%%
clear
FDLabFolders

%% Time-stepping for a generic linear system of ODEs
% *Task*: Implement the $\theta$-scheme for a generic linear system of ODEs.
dbtype FDSolveLinearODEs

%%
% Configure a small test problem whose solution exhibits mild growth
clear
M = [1 1 0; 0 1 0; 0 0 1]; % non-singular mass matrix
A = -[1 1 1; 0 1 2; 0 -2 1]; % RHS system matrix 
n = length(A); % system dimension
y0 = (1 : n)'; % initial states
numSteps = 20; % number of time steps
tGrid = linspace(0, 3, numSteps + 1)'; % grid of time points
[exactSolution, ODESource] = DataLinearODE(M, A, tGrid(1), y0);

%%
% Solution with our new code
theta = 2/3; 
u = FDSolveLinearODEs(M, A, ODESource, tGrid, y0, theta); % solve IVP
%%
% Reference solution with |ode45|
options = odeset('Mass', M);
[~, uref] = ode45(@(t, u) A*u + ODESource(t), tGrid, y0, options); 
%%
% Visual comparison of solutions
clf
set(gca, 'FontSize', FontSize)
plot( ...
    tGrid, u, 'o', ...
    tGrid, uref, '.', ...
    tGrid, FDEvaluateRows(exactSolution, tGrid), '-', ...
    'MarkerSize', 15)
xlabel('t')
ylabel('u(t)')
%%
assert(norm(Compare(u, uref), inf) < 1e-12, 'Mismatch exceeds 1%')
