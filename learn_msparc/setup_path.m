if ispc
    addpath(join([pwd, '\src']))  % windows
else % isunix
    addpath(join([pwd, '/src']))
end

% Then call msparc(prefix) where prefix is input file prefix (string)
% Example:
% msparc('SiH4_quick')