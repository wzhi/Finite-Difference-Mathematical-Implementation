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

