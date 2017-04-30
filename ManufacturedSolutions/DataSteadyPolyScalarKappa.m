function [exact, src, bc] = DataSteadyPolyScalarKappa(kappa, rho)
% Polynomial data for steady PDE.
% [U,SRC,BC]=DATASTEADYPOLYSCALAR(KAPPA,RHO) returns the solution U 
% to the boundary value problem 
%
%  DIV(KAPPA*GRAD(U)) + RHO*U = SRC(T,X1,X2)  on the interior
%        [N1, N2]*GRAD(U) = BC(T,X1,X2,N1,N2) on the boundary
%
% * KAPPA must be a scalar diffision coefficient
% * RHO is a scalar reaction coefficient
%
% The solution U(X1,X2) is a bi-quadratic functon of spatial variables
% (X1, X2), which may be reproduced exactly with 3-point differencing.
%
% See also DATASTEADYNONPOLY, DATAUNSTEADYPOLY.

if ~license('test', 'symbolic_toolbox') || verLessThan('symbolic', '4.0.0')
  % Use previously generated functions
  assert(isequal(kappa, 5))
  assert(isequal(rho, 1))
  % NB: STR2FUNC unavailable before R2009a
  exact = @DataSteadyPolyScalarKappaExactSolution;
  src = @DataSteadyPolyScalarKappaDomainSource;
  bc = @DataSteadyPolyScalarKappaNeumannSource;
  return
end

% If kappa is scalar, convert to multiple of 2x2 identity
kappa = kappa*eye(2); 

% Manufactured solution: An arbitrary biquadratic function
x1 = sym('x1');
x2 = sym('x2');
u = 1 + 2*x1 + 3*x2; % linear part
u = u + x1^2 - 2*x2^2 + 3*x1*x2; % quadratic part

% Generalized Laplacian
gradU = jacobian(u, [x1 x2]).'; % gradient
kappaGradU = kappa*gradU(:);
diffusionTerm = diff(kappaGradU(1), x1) + diff(kappaGradU(2), x2);
reactionTerm = rho*u;
src = diffusionTerm + reactionTerm;

% Normal derivative on the boundary
n1 = sym('n1');
n2 = sym('n2');
bc = n1*gradU(1) + n2*gradU(2);

% Generate regular MATLAB functions
vars = [x1 x2];
exact = FDDataFunction(u, vars, mfilename, 'ExactSolution');
src = FDDataFunction(src, vars, mfilename, 'DomainSource');
bc = FDDataFunction(bc, [vars n1 n2], mfilename, 'NeumannSource');
