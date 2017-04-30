function M = FDMassMatrix(Grid, Dirichlet)
%FDMASSMATRIX Mass matrix for semi-discretized IBVP.
% See also SolveLinearODEs, ODESET.

narginchk(2, 2)
[ny, nx] = size(Grid.Indices);
n = numel(Grid.Indices);
M=speye(n);
for i=1:nx
    M(Grid.Indices(1,i),:)=zeros(1,n);
    M(Grid.Indices(ny,i),:)=zeros(1,n);
end
for i=1:ny
    M(Grid.Indices(i,1),:)=zeros(1,n);
    M(Grid.Indices(i,nx),:)=zeros(1,n);    
end
for i=1:numel(Dirichlet.Indices)
    M(Dirichlet.Indices(i),:)=zeros(1,n);
end