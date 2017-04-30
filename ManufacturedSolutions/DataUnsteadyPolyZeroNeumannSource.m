function expr = DataUnsteadyPolyZeroNeumannSource(t,x1,x2,n1,n2)
%DATAUNSTEADYPOLYZERONEUMANNSOURCE
%    EXPR = DATAUNSTEADYPOLYZERONEUMANNSOURCE(T,X1,X2,N1,N2)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    29-Jul-2016 13:56:10

t2 = cos(t);
t3 = x1-3.0;
t4 = x2-3.0;
expr = n2.*(t2.*t3.*t4.*x1+t2.*t3.*x1.*x2)+n1.*(t2.*t3.*t4.*x2+t2.*t4.*x1.*x2);
