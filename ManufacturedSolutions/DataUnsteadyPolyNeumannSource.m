function expr = DataUnsteadyPolyNeumannSource(t,x1,x2,n1,n2)
%DATAUNSTEADYPOLYNEUMANNSOURCE
%    EXPR = DATAUNSTEADYPOLYNEUMANNSOURCE(T,X1,X2,N1,N2)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    06-Aug-2016 22:00:28

expr = n1.*(x1.*2.0+x2+2.0)+n2.*(x1+x2.*2.0+3.0);
