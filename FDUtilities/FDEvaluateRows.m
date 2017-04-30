function result = FDEvaluateRows(fun, t)
%FDEVALUATEROWS 
% FT=FDEVALUATE(FUN,T) evaluates FUN for each element T(K) of T, 
% storing the results into successive rows of FT. Row I of FT 
% has one column for each element of FUN(GRID(I)).  

narginchk(2, 2)
assert(isnumeric(t))

m = numel(t);
assert(0 < m)
result = flatten(fun(t(1)));
result(m, end) = 0; % pre-allocate storage
for k = 2 : m
  result(k, :) = flatten(fun(t(k)));
end

function x = flatten(x)
x = x(:).';
