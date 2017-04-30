function Y = FDSolveLinearODEs(M, A, f, t, y0, theta, Solver)
%FDSOLVELINEARODES  Solve linear system of ODEs.
% YY=FDSOLVELINEARODES(M,A,F,T,Y0,THETA) returns a numerical
% solution to the initial value problem
%    M*Y'(t) = A*Y(t) + F(t)
%      Y(T(1)) = Y0
%        T(1) <= t <= T(end),
% where
%    M is a constant mass matrix, or [] if no mass matrix is required
%    A is a constant stiffness matrix
%    F is a given function of time
%    Y0 is a given vector containing the initial state
%    [T(1), T(end)] is the time interval of interest
%    YY(k,:) is our numerical approximation of Y(T(k)).
%
% FDSOLVELINEARODES(...,SOLVER) employs a custom algorithm for
% solving the system of linear equations that results at each time step.
%    SOLVER(AA,BB,X0) should compute X such that AA*X=BB
%       where AA is an arbitrary matrix,
%             BB is an arbitrary vector,
%         and X0 is an "initial guess" of X.
%
% Time step sizes are given by DIFF(T)
%  i.e. the size of the K'th step is T(K+1)-T(K).
%
% THETA dictates the time-stepping scheme:
%    THETA   Method           Order of Accuracy
%      0.0   Explicit Euler   1st order
%      1/2   Crank-Nicolson   2nd order
%      2/3   Galerkin Method  1st order
%      1.0   Implicit Euler   1st order
%
% See also ODE45, ODE15S.

narginchk(6, 7)

if nargin < 7 || isempty(Solver)
    % Default solver for linear algebraic equations 
    Solver = @(AA, bb, x0) AA\bb;
end

Y = FDSolveLinearODEsImplementation(M, A, f, t, y0, theta, Solver);
