function MakeAssignment1Report
%MAKEASSIGNMENT1REPORT Joins lab scripts and publishes PDF report.

target = 'Assignment1Report';
files = dir('Task*.m');
format = 'html';

cat = 'type';
if ~strcmp(computer, {'PCWIN', 'PCWIN64'})
    cat = 'cat'; % 'GLNXA64' | 'MACI64'
end

fprintf('Detected %s platform...\n', computer)
fprintf('Using ''%s'' to join files...\n', cat)

addSuffix = @(before, after) strcat(before, '.', after);
fullTarget = addSuffix(target, 'm');
try
    command = sprintf('%s %s > %s', ...
        cat, strjoin({files.name}, ' '), fullTarget);
    fprintf('Executing ''%s''...\n', command)
    assert(system(command) == 0, 'Command failed: "%s", command')
catch e
    fprintf('Error: %s\n', e.message)
    disp('Please check your files.')
    disp('In a real emergency, upload your zip file and email Jon.')
    return
end

fprintf('Opening %s...\n', fullTarget)
edit(fullTarget)

publish(target, format)

report = fullfile('.', 'html', addSuffix(target, format));
fprintf('Report saved to %s\n', report)
try
    open(report)
catch e
    disp(e.message)
end
