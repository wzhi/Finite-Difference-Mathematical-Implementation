function [U, A, b] = FDSolve( ...
    Grid, kappa, rho, domainSource, NeumannSource, Dirichlet)
%FDSOLVE Finite difference solution of reaction-diffusion equation.
% See also FDSYSTEMMATRIX, FDSYSTEMVECTOR.

timer = StartTimer('system matrix');
A = FDSystemMatrix(Grid, kappa, rho, Dirichlet);
StopTimer(timer)

timer = StartTimer('system RHS');
b = FDSystemVector(Grid, domainSource, NeumannSource, Dirichlet);
StopTimer(timer)

timer = StartTimer('solve');
U = reshape(A\b, size(Grid.Indices));
StopTimer(timer)
