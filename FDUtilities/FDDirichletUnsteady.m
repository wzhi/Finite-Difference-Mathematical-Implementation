function Dirichlet = FDDirichletUnsteady(grid, fun, varargin)
%FDDIRICHLETUNSTEADY  Unsteady Dirichlet/Type-1 BC for 2-dimensional grid.
% See also FDGRID, FDDIRICHLET.

narginchk(2, inf)

% Concatenate lists of indices in optional arguments
varargin = cellfun(@(x) x(:), varargin, 'UniformOutput', false);
indices = unique(vertcat(varargin{:}));
unsteady = @(t) fun(t, grid.X(indices), grid.Y(indices));

Dirichlet = struct('Indices', indices(:), 'RHSValues', unsteady);
