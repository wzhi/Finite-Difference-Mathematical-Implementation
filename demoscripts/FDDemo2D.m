%% Numerical solution of a Poisson problem on the unit square

%% Overview
% Numerical solution of the Poisson equation with Dirichlet boundary
% condition on the unit square:
%
% $\nabla^2 u = $ domainSource $\quad\leftarrow$ on unit square
%
% $\quad u = $ boundarySource $\quad\leftarrow$ on boundary
%

%% RHS/Source functions
syms x y
u = x^2 + 3*x*y + 2*x - 2*y^2 + 3*y;
exactSolution = matlabFunction(u); %#ok<*NOPTS>
domainSource = matlabFunction(diff(u, x, 2) + diff(u, y, 2), 'vars', [x, y]);
boundarySource = exactSolution; % because of Dirichlet BCs!

%% Grid
nx = 4; ny = 3; % number of grid points
x = linspace(0, 1, nx); % "ticks" on grid axes
y = linspace(0, 1, ny); %
hx = (x(end) - x(1))/(nx - 1); % grid spacing
hy = (y(end) - y(1))/(ny - 1);

%%
[X, Y] = meshgrid(x, y); % grid-point coordinate pairs
grid2ind = reshape(1 : nx*ny, size(X)); % grid-point indices
%%
disp([grid2ind, X, Y]) % concatenate arrays merely to save page space

%%
clf
plot(X, Y, ...
    'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r'), hold on
text(X(:) + 0.01, Y(:), num2str(grid2ind(:)), 'FontSize', 15)
set(gca, 'XTick', x, 'YTick', y)
xlabel('x'), ylabel('y'), axis equal, axis tight

%% Local matrix of finite difference coefficients
LaplacianMolecule = [
    0               1/hy^2             0;
    1/hx^2  -2*(1/hx^2 + 1/hy^2)  1/hx^2;
    0               1/hy^2             0;
    ]

%% Assemble discrete equations
A = eye(numel(grid2ind)); % global matrix of coefficients
b = zeros(numel(grid2ind), 1); % global left-hand side matrix
for i = 2 : size(grid2ind, 1) - 1
    for j = 2 : size(grid2ind, 2) - 1
        %%
        fprintf(' ---- grid point %d, %d ----\n', i, j)
        centre = grid2ind(i, j);
        stencil = grid2ind(i-1 : i+1, j-1 : j+1)
        A(centre, stencil) = LaplacianMolecule(:);
        b(centre) = domainSource(X(centre), Y(centre));        
    end
end

%% RHS entries of BC equations
b(grid2ind(:, [1 end])) = boundarySource(X(:, [1 end]), Y(:, [1 end]));
b(grid2ind([1 end], :)) = boundarySource(X([1 end], :), Y([1 end], :));
%%
labels = 1 : numel(X); % row/column labels
disp([ nan, labels, nan; labels', A, b]) % concatenate to save page space

%% Solve and plot
U = reshape(A\b, size(X));
%%
clf, hold on
[XX, YY] = meshgrid(linspace(x(1), x(end), 10), linspace(y(1), y(end), 10));
surf(XX, YY, exactSolution(XX, YY), 'FaceAlpha', 0.5);
plot3(X, Y, double(U), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
xlabel('x'), ylabel('y'), view(3), grid on

%% Verify solution
assert(norm(exactSolution(X, Y) - U, inf) < 1e-12)
