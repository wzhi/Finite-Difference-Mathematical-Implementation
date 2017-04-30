function [y, f] = DataLinearODE(M, A, t0, y0)
%DATALINEARODE  Defines system of linear ODEs with time-linear solution.
% [U,F]=DATAODEQUADRATIC(M,A,T0,Y0) returns the exact solution U(T) to
% the initial value problem:
%      M*U'(T) = A*U + F(T)
%        U(T0) = U0

n = length(A);
v = (-1).^mod(1 : n, 2)'; % [+1; -1; +1; ... ]
y = @(t) y0 + v*(t - t0); 
dydt = @(~) v;
f = @(t) M*dydt(t) - A*y(t);
