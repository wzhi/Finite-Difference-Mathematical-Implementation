function expr = DataSteadyPolyExactSolution(x1,x2)
%DATASTEADYPOLYEXACTSOLUTION
%    EXPR = DATASTEADYPOLYEXACTSOLUTION(X1,X2)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    10-Aug-2016 23:06:59

expr = x1.*2.0+x2.*3.0+x1.*x2.*3.0+x1.^2-x2.^2.*2.0+1.0;
