function [exact, src, bc] = DataUnsteadyNonPoly(kappa, rho)
%DATAUNSTEADYPOLY  Non-polynomial data for unsteady PDE.
% [U,SRC,BC]=DATAUNSTEADYNONPOLY(KAPPA,RHO) returns the solution U to the
% initial boundary value problem
%
%  DIV(KAPPA*GRAD(U)) + RHO*U - D(U)/DT = SRC(T,X1,X2)  on the interior
%                  [N1, N2]*GRAD(U) = BC(T,X1,X2,N1,N2) on the boundary
%
% The solution U(T,X1,X2) cannot be reproduced exactly with central differencing. 
%
% See also DATASTEADYPOLY, DATASTEADYNONPOLY.

if ~license('test', 'symbolic_toolbox') || verLessThan('symbolic', '4.0.0')
  % Use previously generated functions
  assert(isequal(kappa, 5))
  assert(isequal(rho, 1))
  % NB: STR2FUNC unavailable before R2009a
  exact = @DataUnsteadyNonPolyExactSolution;
  src = @DataUnsteadyNonPolyDomainSource;
  bc = @DataUnsteadyNonPolyNeumannSource;
  return
end

% If kappa is scalar, convert to multiple of 2x2 identity
kappa = kappa*eye(2);

% Manufactured solution
t = sym('t');
x1 = sym('x1');
x2 = sym('x2');
u = sin(t) + sin(x1) + cos(x2) - sin(x2)*cos(x1);

% Generalized Laplacian
du = jacobian(u, [x1 x2]).';
adu = kappa*du(:);
lapu = diff(adu(1), x1) + diff(adu(2), x2);
src = lapu + rho*u - diff(u, t);

% Normal derivative on the boundary
n1 = sym('n1');
n2 = sym('n2');
bc = n1*du(1) + n2*du(2);

% Generate regular MATLAB functions
vars = [t x1 x2];
exact = FDDataFunction(u, vars, mfilename, 'ExactSolution');
src = FDDataFunction(src, vars, mfilename, 'DomainSource');
bc = FDDataFunction(bc, [vars n1 n2], mfilename, 'NeumannSource');

