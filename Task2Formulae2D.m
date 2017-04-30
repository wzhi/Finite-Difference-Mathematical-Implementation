%% Applying FD formulae in 2D

%%
clear
FDLabFolders
tol = 1e-12;

%% Define target polynomial and its derivatives
FDTask2Data % defines polynomial and its derivatives
p2Analytical = p2([xstar; ystar]);

%% Evaluate 1D finite difference coefficients
cx0 = FDMolecule1D(0, xstar, x); % <== TODO ==
cx1 = FDMolecule1D(1, xstar, x); % <== TODO ==
cx2 = FDMolecule1D(2, xstar, x); % <== TODO ==

cy0 = FDMolecule1D(0, ystar, y); % <== TODO ==
cy1 = FDMolecule1D(1, ystar, y); % <== TODO ==
cy2 = FDMolecule1D(2, ystar, y); % <== TODO ==

%% Apply coefficients to function values on 2D grid
d2pdx2 = cy0(:)'*U*cx2(:);
d2pdy2 = cy2(:)'*U*cx0(:);
d2pdxdy = cy1(:)'*U*cx1(:);

%% Test second derivatives
p2Numerical = [d2pdx2, d2pdxdy; d2pdxdy, d2pdy2];
assert(norm(Compare(p2Analytical, p2Numerical), inf) < tol)

%% 2D FD molecule for mixed second derivative
% Here, we compute a 3-by-3 molecule (3-point in each direction), which is
% comparable with |LaplacianMolecule| in |FDDemo2D.m|.
B = cx1*cy1';
disp(B);
%%
d2pdxdy = B(:)'*U(:);
assert(Compare(p2Analytical(1, 2), d2pdxdy) < tol)
