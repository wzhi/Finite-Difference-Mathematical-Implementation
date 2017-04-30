function func = FDDataFunction(expr, vars, mfile, description, path)
%FDDATAFUNCTION  Creates M-file from symbolic expression.
% Stores the resulting M-file in the Data directory.
if nargin < 5 || isempty(path)
  path = pwd;
end
subfolder = 'ManufacturedSolutions';
destination = fullfile(path, subfolder, [mfile description '.m']);
func = matlabFunction(expr, 'vars', vars, 'file', destination);
