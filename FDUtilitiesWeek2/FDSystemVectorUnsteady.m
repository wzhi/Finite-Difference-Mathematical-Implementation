function b = FDSystemVectorUnsteady(t, ...
    Grid, domainSource, NeumannSource, Dirichlet)
%FDSYSTEMVECTORUNSTEADY Adapts FDSystemVector for time-varying data.
% See documentation for FDSYSTEMVECTOR and FDSYSTEMMATRIX.

% Modify interface of time-varying functions to appear "steady"
domainSourceSteady = @(x, y) domainSource(t, x, y);
NeumannSourceSteady = @(x, y, nx, ny) NeumannSource(t, x, y, nx, ny);
if ~isnumeric(Dirichlet.RHSValues)
    % Dirichlet data may (optionally) be time-varying
    Dirichlet.RHSValues = Dirichlet.RHSValues(t);
end

% Delegate to the steady version
b = FDSystemVector(Grid, domainSourceSteady, NeumannSourceSteady, Dirichlet);
