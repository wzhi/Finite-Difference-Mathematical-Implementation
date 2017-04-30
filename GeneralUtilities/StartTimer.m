function timer = StartTimer(label, fid) 
%STARTTIMER Start a stopwatch with given label.
% See also STOPTIMER, TIC, TOC.
narginchk(1, 2)
if nargin < 2 || isempty(fid)
    fid = 1; % defaults to stdout
end
assert(ischar(label))
fprintf('Starting %s... \n', label)
timer = struct( ...
    'FileID', fid, ...
    'Label', label, ...
    'StartTime', builtin('tic'));
