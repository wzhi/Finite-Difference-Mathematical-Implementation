%% Numerical solution of a reaction-diffusion problem on the unit interval

%% Overview
% We are aiming to find a numerical solution to a
% 1D reaction-diffusion equation on the unit interval. 
%
% $\quad u + \nabla^2 u = $ domainSource $\quad\leftarrow$ on unit interval
%
% $\quad n \cdot \nabla u = $ NeumannSource $\quad\leftarrow$ at west/left boundary
%
% $\quad u = $ DirichletSource $\quad\leftarrow$ at east/right boundary

%% RHS/Source functions
syms x n
u = 1 - 3*x + 4*x.^2; % prescribed solution
exactSolution = matlabFunction(u); %#ok<*NOPTS>
domainSource = matlabFunction(diff(u, x, 2) + u);
NeumannSource = matlabFunction(n*diff(u, x), 'vars', [x, n]);
DirichletSource = exactSolution; % because of Dirichlet BCs!

%% Grid
% The expressions here are best motivated by a comparison with the 2D case.
nx = 5; % number of grid points
x = linspace(0, 1, nx) % grid points on x-axis...
hx = (x(end) - x(1))/(nx - 1); % grid spacing
%%
X = x; % ... comprise the grid itself (compare with 2D case)
grid2ind = 1 : numel(X); % grid-point indices (trivial in 1D)
disp([grid2ind; X]) % concatenate arrays merely to save page space

%% Local matrix of finite difference coefficients
ReactionMolecule = [0, 1, 0];
DiffusionMolecule = [1, -2, 1]/hx^2

%% Assemble discrete equations
A = zeros(numel(grid2ind)); % global matrix of coefficients
b = zeros(numel(grid2ind), 1); % global left-hand side matrix
for i = grid2ind(2 : end-1)
    centre = grid2ind(i);
    stencil = grid2ind(i-1 : i+1);
    fprintf('centre: %d --> stencil: %s\n', centre, mat2str(stencil))
    A(centre, stencil) = ReactionMolecule(:) + DiffusionMolecule(:);
    b(centre) = domainSource(X(centre));
end

%% Neumann/flux condition at left boundary
Molecule = [-3, 4, -1]/(2*hx);
OutwardNormal = -1;
A(grid2ind(1), grid2ind(1 : 3)) = OutwardNormal*Molecule;
b(grid2ind(1)) = NeumannSource(X(1), OutwardNormal);

%% Dirichlet condition at right boundary
A(grid2ind(end), grid2ind(end)) = 1;
b(grid2ind(end)) = DirichletSource(X(end));
%%
labels = 1 : numel(X); % row/column labels
disp([ nan, labels, nan; labels', A, b]) % concatenate to save page space

%% Solve and plot
U = reshape(A\b, size(X)); % compare with 2D case
xx = linspace(x(1), x(end), 100);
plot(X, U, 'o', xx, exactSolution(xx), '-')
set(legend('FD', 'exact', 'Location', 'Best'), 'box', 'off')

%% Verify solution
assert(norm(U - exactSolution(X)) < 1e-12)
