function ftxy = FDEvaluate(fun, grid, tgrid)
%FDEVALUATE Evaluate function on spatio-temporal grid.
%
% FDEVALUATE(FUN,GRID) evaluates FUN(GRID.X,GRID.Y) on a spatial GRID
% returned by FDGRID. The result is the same shape as GRID.IND.
% Note that FUN Must be vectorized.
%
% FDEVALUATE(FUN,GRID,TGRID) evaluates FUN(T,X,Y) on a spatio-temporal 
% grid defined by numeric arrays TGRID, GRID.X, and GRID.Y. 
% The result has one row for each temporal point and one column for each
% spatial point, consistent with MATLAB's ODE Suite.  
%
narginchk(2, 3)
assert(isstruct(grid))
if nargin == 2 
    % Purely spatial grid
    ftxy = fun(grid.X, grid.Y);
else
    % Spatio-temporal grid
    ftxy = FDEvaluateRows(@(t) fun(t, grid.X, grid.Y), tgrid);
end
