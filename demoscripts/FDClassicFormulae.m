%% Classical finite difference formulae

%%
clear
h = 0.1; % grid spacing
x = h*(1 : 5); % uniformly spaced grid
Display = @(d, k, n, c) ...
    fprintf('for u^(%d) @ x(%d) of x(1 : %d): %s/(%g h^%d)\n', ...
    d, k, n, mat2str(FDMolecule1D(d, x(k), x(1 : n))*c*h^d, 2), c, d);

%% Table of Standard 1st-order Formulae 
% See Section 3.1 of the Notes
Display(1, 1, 2, 1)
Display(1, 2, 2, 1)
Display(1, 2, 3, 2)
Display(1, 1, 3, 2)
Display(1, 3, 3, 2)
Display(1, 2, 4, 6)
Display(1, 3, 4, 6)
Display(1, 3, 5, 24)

%% Table of 2nd-order formulae 
% See Section 3.2 of the Notes
Display(2, 2, 3, 1)
Display(2, 1, 3, 1)
Display(2, 3, 3, 1)
Display(2, 3, 5, 12)
