function b = FDSystemVector(Grid, domainSource, NeumannSource, Dirichlet)
%FDSYSTEMVECTOR  Right-hand side vector corresponding to FDSYSTEMMATRIX.
% B=FDSYSTEMVECTOR(GRID,SRC,BC,DIRICHLET) generates the global right-hand
% side vector for the 3-point finite difference discretization of the
% problem
%       GRAD(AD*GRAD(U)) + AR*U = SRC(X,Y) on the interior of (X,Y)
%                     N*GRAD(U) = BC(X,Y) on the boundary of (X,Y)
%          U(DIRICHLET.Indices) = DIRICHLET.RHSValues
%
% See also FDGRID, FDDIRICHLET, FDSYSTEMMATRIX.

b = zeros(numel(Grid.Indices), 1); % pre-allocate result


%
% TODO, Steps 2 to 4: Assign to rows corresponding to all other grid points.
%       Step 1 is already finished, below.
%

% Dirichlet points
[ny, nx] = size(Grid.Indices);
for i = 1 : ny
    b(Grid.Indices(i,1))=NeumannSource(Grid.X(i,1), Grid.Y(i,1), -1, 0);
    b(Grid.Indices(i,nx))=NeumannSource(Grid.X(i,nx),Grid.Y(i,nx),1, 0);
end
for j = 2 : nx - 1
    b(Grid.Indices(1,j))=NeumannSource(Grid.X(1,j),Grid.Y(1,j),0,-1);
    b(Grid.Indices(ny,j))=NeumannSource(Grid.X(ny,j),Grid.Y(ny,j),0, 1);
end

for i = 2 : ny - 1
    for j = 2 : nx - 1
        b(Grid.Indices(i,j))=domainSource(Grid.X(i,j),Grid.Y(i,j));
    end
end
b(Dirichlet.Indices) = Dirichlet.RHSValues;

