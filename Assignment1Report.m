%% Computing and applying FD formulae in 1D

%%
clear
FDLabFolders
tol = 1e-14; % (there is no truncation error)

%% Task 1(a)
ToDo('Study- and implement FDMolecule1D')
dbtype FDMolecule1D % <== TODO ==

%% Task 1(b)
% Comparison with |ppval| and |ppder|
x = [0 1 3 5 6]; % grid points
xstar = 4; % evaluation point
%%
% You may well use more than one line for each calculation!
% Analytical:
analytical=polyval(polyder(polyder(1:5)),xstar);
%F'(XSTAR) = DOT(F(X), C)
c=FDMolecule1D(2, xstar, x)
numerical = dot(polyval((1:5),x),c); % <== TODO ==
%%
fprintf('analytical: %g\n numerical: %g\n', analytical, numerical)
assert(Compare(analytical, numerical) < tol)
assert(Compare(analytical, 246) < tol) % This is the correct value!
%% Applying FD formulae in 2D

%%
clear
FDLabFolders
tol = 1e-12;

%% Define target polynomial and its derivatives
FDTask2Data % defines polynomial and its derivatives
p2Analytical = p2([xstar; ystar]);

%% Evaluate 1D finite difference coefficients
cx0 = FDMolecule1D(0, xstar, x); % <== TODO ==
cx1 = FDMolecule1D(1, xstar, x); % <== TODO ==
cx2 = FDMolecule1D(2, xstar, x); % <== TODO ==

cy0 = FDMolecule1D(0, ystar, y); % <== TODO ==
cy1 = FDMolecule1D(1, ystar, y); % <== TODO ==
cy2 = FDMolecule1D(2, ystar, y); % <== TODO ==

%% Apply coefficients to function values on 2D grid
d2pdx2 = cy0(:)'*U*cx2(:);
d2pdy2 = cy2(:)'*U*cx0(:);
d2pdxdy = cy1(:)'*U*cx1(:);

%% Test second derivatives
p2Numerical = [d2pdx2, d2pdxdy; d2pdxdy, d2pdy2];
assert(norm(Compare(p2Analytical, p2Numerical), inf) < tol)

%% 2D FD molecule for mixed second derivative
% Here, we compute a 3-by-3 molecule (3-point in each direction), which is
% comparable with |LaplacianMolecule| in |FDDemo2D.m|.
B = cx1*cy1';
disp(B);
%%
d2pdxdy = B(:)'*U(:);
assert(Compare(p2Analytical(1, 2), d2pdxdy) < tol)
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
%% Mass matrix for semi-discretized IBVPs

%%
clear
FDLabFolders

%%
dbtype FDMassMatrix
%%
xyGrid = FDGrid(0 : 3);
Dirichlet = FDDirichlet( ...
    xyGrid, @(t, x, y) zeros(size(x)), ...
    xyGrid.Indices([1 end], :)); % points on north & south boundaries, say
%%
% "Does this output meet with your expectation?"
disp(FDMassMatrix(xyGrid, Dirichlet))
%% Harbour Analysis

%%
clear
FDLabFolders

%%
secondsPerDay = 60*60*24;
rho = 0; % reaction coefficient
kappa = 7; % diffusion coefficient
%% our functions for domain/neumann source
myDomainSource=@(t,x,y) zeros(size(x));
myNeumannSource=@(t,x,y,nx,ny) zeros(size(x)); %All zeros anyway.

%%
% Spatial grid
x = 0:50:1000;
y = 0:50:3000;
xyGrid = FDGrid(x, y); % xyGrid point coordinates

opening=unique(interp2(xyGrid.X,xyGrid.Y,xyGrid.Indices,400:600,3000,'nearest'));%get indices of opening
rhsOpening=zeros(1,numel(opening));
xTargets = [1000, 1000];
yTargets = [2300, 1000];
%%
% Indices of target points
targetIndices = interp2( ...
    xyGrid.X, xyGrid.Y, xyGrid.Indices, ...
    xTargets, yTargets, 'nearest');
indices=[targetIndices opening];
rhs=[1.3 1.3 rhsOpening];
Dirichlet = struct('Indices', indices(1,:), 'RHSValues', rhs(1,:));

%%
% Temporal grid
timeSpan = [0 30*secondsPerDay]; % time interval
numSteps = 300; % number of time steps
tGrid = linspace(timeSpan(1), timeSpan(end), numSteps + 1)';

%%
% Initial state
uInitial = zeros(size(xyGrid.Indices));
uInitial(targetIndices)=1.3; %initial conditions of pipes set at 1.3
%%
% Mass matrix and stiffness matrix
M = FDMassMatrix(xyGrid, Dirichlet);
A = FDSystemMatrix(xyGrid, kappa, rho, Dirichlet);

%%
% Modify interface for use with generic time-stepping routines
f = @(t) -FDSystemVectorUnsteady( ...
    t, xyGrid, myDomainSource, myNeumannSource, Dirichlet);

%%
% Solution with our new code
theta = 0.5;
U = FDSolveLinearODEs(M, A, f, tGrid, uInitial, theta); % solve IVP
hold on;

xSamples = [100, 500, 600];
ySamples = [1500, 2500,100];


% Indices of target points
%interpolate for samples
samples = interp2( ...
    xyGrid.X, xyGrid.Y, xyGrid.Indices, ...
    xSamples, ySamples, 'nearest');
marker = {'MarkerSize', 3};
plot(xyGrid.X, xyGrid.Y, 'k.', marker{:});
plot(xSamples,ySamples,'ro', marker{:});
plot(400:600,3000,'bs', marker{:});
plot(xTargets,yTargets,'bs', marker{:});
axis equal, axis tight
hold off;
%%
plot(0:0.1:30,U(:,samples(1)),0:0.1:30,U(:,samples(2)),0:0.1:30,U(:,samples(3)));
%making the surf:
%%
z=zeros(61,21);
for i=1:61
    for j=1:21
        z(i,j)=U(end,xyGrid.Indices(i,j));
    end
end

surf(x,y,z);
%%

%% Reports

%% Part 6(a)
%
% 1.)Formulae for |AA| and |bb|
%
% $AA=\frac{1}{\Delta t}M-\theta A$
%
% $bb=(\frac{1}{\Delta t}M+\overline{\theta} A)U_n + f(\overline{\theta}t_n+\theta t_{n+1})$
%
% 2.)\ is currently used. This will result in LU factorisation and
% substitution to arive at a solution directly. The upper and lower matrices are not stored. 
% We could use an iterative
% method, by using the X0 inputted. This may result in less computations
% required to arrive at a numerical solution. Especially, if we select
% good initial conditions.
%

%% Part 6(b)
%
% 1.) 
%
%           -We have approximated a 3d problem to 2d.  
%
%           -The boundary conditions of the open ocean is an unaccurate
%           approximation, as the poluutants there will not all disperse
%           immediately. It is unreasonable hold it constant at 0.
%
%           -We have assumed that the concentration of the polutant at
%           where the pipes are located have a constant concentration. In
%           reality, the concentration at these points may not be constant,
%           as pollutants are diluted instantly when in water.
%
%           -The Neumann conditions may also not be realistic, as it may be
%           a function, and not held constant.
%
%           -The kappa term, diffusion, may not be a scalar and instead be
%           a tensor.
%
% 2.)At all three of the sample points, the concentration is below
% 1.0kgm^-1. Hence, waste concentration requirements are met. All three
% samples have concentrations that have stabilised by 20 days. Taking the
% sample at 20 days is sensible.
%
% 3.) We could compute multiple times with decreasing step size, and
% compare results. If the difference is within a decided tolerance, then we can
% terminate.
%
