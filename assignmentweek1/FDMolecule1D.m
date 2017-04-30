function c = FDMolecule1D(d, xstar, x)
%FDMOLECULE1D Finite difference moledule on a 1D grid.
% FDMOLECULE1D(D,XSTAR,X) generates the NUMEL(X)-point finite 
% difference coefficients for approximating the D'th derivative at XSTAR of
% a function sampled at X 
%      e.g. if C = FDMOLEDULE1D(1, XSTAR, X),
%         then F'(XSTAR) = DOT(F(X), C) + "truncation error".
% See also FDSYSTEMMATRIX, FDSYSTEMVECTOR.

narginchk(3, 3)
assert(d < numel(x))
NML=numel(x);
A=ones(NML,NML);
b=zeros(NML,1);
b(d+1,1)=factorial(d);
for i=2:NML
    for j=1:NML
        A(i,j)=(x(j)-xstar)^(i-1);
        
    end
end
c=A\b;
