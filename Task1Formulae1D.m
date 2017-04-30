%% Computing and applying FD formulae in 1D

%%
clear
FDLabFolders
tol = 1e-14; % (there is no truncation error)

%% Task 1(a)
ToDo('Study- and implement FDMolecule1D')
dbtype FDMolecule1D % <== TODO ==

%% Task 1(b)
% Comparison with |ppval| and |ppder|
x = [0 1 3 5 6]; % grid points
xstar = 4; % evaluation point
%%
% You may well use more than one line for each calculation!
% Analytical:
analytical=polyval(polyder(polyder(1:5)),xstar);
%F'(XSTAR) = DOT(F(X), C)
c=FDMolecule1D(2, xstar, x)
numerical = dot(polyval((1:5),x),c); % <== TODO ==
%%
fprintf('analytical: %g\n numerical: %g\n', analytical, numerical)
assert(Compare(analytical, numerical) < tol)
assert(Compare(analytical, 246) < tol) % This is the correct value!
