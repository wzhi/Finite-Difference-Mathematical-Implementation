function varargout = StopTimer(timer)
%STOPTIMER Stops a stopwatch with message.
% See also STARTTIMER, TIC, TOC.
narginchk(1, 1)
elapsed = builtin('toc', timer.StartTime);
if nargout == 1
    varargout{1} = elapsed;
else
    fprintf(timer.FileID, ...
        '     ... %s finished: %1.2g sec\n', timer.Label, elapsed);
end
