function expr = DataSteadyPolyNeumannSource(x1,x2,n1,n2)
%DATASTEADYPOLYNEUMANNSOURCE
%    EXPR = DATASTEADYPOLYNEUMANNSOURCE(X1,X2,N1,N2)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    10-Aug-2016 23:06:59

expr = n1.*(x1.*2.0+x2.*3.0+2.0)+n2.*(x1.*3.0-x2.*4.0+3.0);