clear all; close all;

if ispc
    delim_char = '\';  % windows
else % isunix
    delim_char = '/';
end
% It seems that addpath will also recognize / as path separator in Windows

addpath(join([pwd, delim_char, 'src']));

% For other functions
% These will be executed by msparc function
addpath(join([pwd, delim_char, 'src', delim_char, 'xc']));
addpath(join([pwd, delim_char, 'src', delim_char, 'xc', delim_char, 'exx']));
addpath(join([pwd, delim_char, 'src', delim_char, 'xc', delim_char, 'mgga']));
addpath(join([pwd, delim_char, 'src', delim_char, 'xc', delim_char, 'vdW']));

% Then call msparc(prefix) where prefix is input file prefix (string)
% Example:
% msparc('SiH4_quick')

