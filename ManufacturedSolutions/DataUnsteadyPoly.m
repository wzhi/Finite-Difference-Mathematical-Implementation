function [exact, src, bc] = DataUnsteadyPoly(kappa, rho)
%DATAUNSTEADYPOLY  Polynomial data for unsteady PDE.
% [U,SRC,BC]=DATAUNSTEADYPOLY(KAPPA,RHO) returns the solution U to the
% initial boundary value problem
%
%  DIV(KAPPA*GRAD(U)) + RHO*U - D(U)/DT = SRC(T,X1,X2)  on the interior
%                  [N1, N2]*GRAD(U) = BC(T,X1,X2,N1,N2) on the boundary
%
% The solution U(T,X1,X2) is a bi-quadratic functon of spatial variables
% (X1, X2), which may be reproduced exactly with central differencing.
%
% See also DATASTEADYPOLY, DATASTEADYNONPOLY.

if ~license('test', 'symbolic_toolbox') || verLessThan('symbolic', '4.0.0')
  % Use previously generated functions
  assert(isequal(kappa, 5))
  assert(isequal(rho, 1))
  % NB: STR2FUNC unavailable before R2009a
  exact = @DataUnsteadyPolyExactSolution;
  src = @DataUnsteadyPolyDomainSource;
  bc = @DataUnsteadyPolyNeumannSource;
  return
end

% If kappa is scalar, convert to multiple of 2x2 identity
kappa = kappa*eye(2);

% Manufactured solution
t = sym('t');
x1 = sym('x1');
x2 = sym('x2');
u = cos(t) + 2*x1 + 3*x2 + x1*x2 + x1^2 + x2^2;

% Generalized Laplacian
dudt = diff(u, t);
gradu = jacobian(u, [x1 x2]).';
Adu = kappa*gradu(:);
lapu = diff(Adu(1), x1) + diff(Adu(2), x2);
f = -dudt + (lapu + rho*u);

% Normal derivative on the boundary
n1 = sym('n1');
n2 = sym('n2');
g = n1*gradu(1) + n2*gradu(2);

% Generate regular MATLAB functions
vars = [t x1 x2];
exact = FDDataFunction(u, vars, mfilename, 'ExactSolution');
src = FDDataFunction(f, vars, mfilename, 'DomainSource');
bc = FDDataFunction(g, [vars n1 n2], mfilename, 'NeumannSource');
