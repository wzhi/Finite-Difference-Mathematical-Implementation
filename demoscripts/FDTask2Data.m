%% This script generates data for Part 2, Week 1

%% Parameters
Q = [1 -3; -2 2]; % symmetric matrix of coefficients of quadratic form 
b = [2 -2]; % coefficients on linear part
c = 3; % constant offset

%% Bi-quadratic polynomial and its derivatives
symm = @(A) 0.5*(A + A'); % matrix symmetrization
p0 = @(x) c + b(:).'*x(:) + x(:).'*Q*x(:); % 0th derivative
p1 = @(x) b(:).' + 2*x(:).'*symm(Q); % 1st derivative
p2 = @(x) 2*symm(Q); % 2nd derivative

%% Finite difference data
% Grid
x = [0 0.5 1.5];
y = [1 2 4];
[X, Y] = meshgrid(x, y);
%%
% Function values on grid
U = arrayfun(@(x, y) p0([x; y]), X, Y);
%% 
% Coordinates of evaluation point
xstar = x(2); 
ystar = y(2);
