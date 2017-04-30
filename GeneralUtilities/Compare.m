function z = Compare(x, y)
%COMPARE  Normalised difference.
% COMPARE(X,Y) compares corresponding components of X and Y. 
% Each component of the result behaves like
% * a relative difference if either argument is large in magnitude, or
% * an absolute difference if both arguments are small in magnitude.
% X and Y must have the same shape, unless one is scalar.
%
% See the ENGSCI 233 notes or the following link for a discussion:
%  https://en.wikipedia.org/wiki/Relative_change_and_difference
%
narginchk(2, 2)
scale = 0.5*(abs(x) + abs(y));
z = abs(x - y)./(1 + scale);
