function expr = DataUnsteadyPolyZeroDomainSource(t,x1,x2)
%DATAUNSTEADYPOLYZERODOMAINSOURCE
%    EXPR = DATAUNSTEADYPOLYZERODOMAINSOURCE(T,X1,X2)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    29-Jul-2016 13:56:10

t2 = cos(t);
t3 = x1-3.0;
t4 = x2-3.0;
expr = t2.*t3.*x1.*1.0e1+t2.*t4.*x2.*1.0e1+t3.*t4.*x1.*x2.*sin(t)+t2.*t3.*t4.*x1.*x2;
