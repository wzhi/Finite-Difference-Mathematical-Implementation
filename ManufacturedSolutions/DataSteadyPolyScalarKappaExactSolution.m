function expr = DataSteadyPolyScalarKappaExactSolution(x1,x2)
%DATASTEADYPOLYSCALARKAPPAEXACTSOLUTION
%    EXPR = DATASTEADYPOLYSCALARKAPPAEXACTSOLUTION(X1,X2)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    06-Aug-2016 01:09:44

expr = x1.*2.0+x2.*3.0+x1.*x2.*3.0+x1.^2-x2.^2.*2.0+1.0;