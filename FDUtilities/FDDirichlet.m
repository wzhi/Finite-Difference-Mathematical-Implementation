function Dirichlet = FDDirichlet(grid, funcOrValues, varargin)
%FDDIRICHLET  Dirichlet/Type-1 BC for 2-dimensional grid.
%
% This function may be used in two ways:
%
% Second argument is numeric
% ==========================
% FDDirichletBC(GRID,VALUES,INDICES)
% conditions of the form
%      U(K) = VALUES(K)
% for each index K in INDICES.
%
% Second argument is a function
% =============================
% FDDirichletBC(GRID,FUN,IND1,IND2,..) encodes conditions of the form
%      U(K) = F(GRID.x(K),GRID.y(K))
% for each index K in [IND1(:); IND2(:); ..].
%
% See also FDGRID.

narginchk(2, inf)

if isnumeric(funcOrValues)
    % Second argument is numeric
    % ==========================
    assert(numel(varargin) == 1)
    indices = varargin{1};
    assert(numel(indices) == numel(funcOrValues))
    values = funcOrValues;
else
    % Secnd argument is a function
    % ============================
    % Concatenate lists of indices in optional arguments
    varargin = cellfun(@(x) x(:), varargin, 'UniformOutput', false);
    indices = unique(vertcat(varargin{:}));
    values = funcOrValues(grid.X(indices), grid.Y(indices));
end

Dirichlet = struct('Indices', indices(:), 'RHSValues', values(:));
