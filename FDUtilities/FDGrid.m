function grid = FDGrid(x, y)
%FDGRID  2-dimensional finite difference grid with given points.
% GRID=FDGRID(X,Y) encodes a two-dimensional grid with points X(:) 
% on the horizontal axis and Y(:) on the vertical axis.
% FDGRID(X) is equivalent to FDGRID(X,X).
% For each grid point (I,J):
% * GRID.X(I,J) = horizontal coordinate
% * GRID.Y(I,J) = vertical coordinate
% * GRID.IND(I,J) = index
% See also FDEvaluate.
narginchk(1, 2)
if nargin < 2 || isempty(y)
    y = x;
end
[x, y] = meshgrid(sort(x), sort(y));
ind = reshape(1 : numel(x), size(x));
grid = struct('X', x, 'Y', y, 'Indices', ind);
