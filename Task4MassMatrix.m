%% Mass matrix for semi-discretized IBVPs

%%
clear
FDLabFolders

%%
dbtype FDMassMatrix
%%
xyGrid = FDGrid(0 : 3);
Dirichlet = FDDirichlet( ...
    xyGrid, @(t, x, y) zeros(size(x)), ...
    xyGrid.Indices([1 end], :)); % points on north & south boundaries, say
%%
% "Does this output meet with your expectation?"
disp(FDMassMatrix(xyGrid, Dirichlet))
