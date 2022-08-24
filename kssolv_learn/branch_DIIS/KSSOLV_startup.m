function KSSOLV_startup()
% KSSOLV_STARTUP  Startup file for KSSOLV
%   MAKE adds paths of the KSSOLV to Matlab.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

file_path = '/home/efefer/WORKS/KSSOLV/kssolv2_DIIS/';

% Folder for all utility functions
addpath([file_path 'util']);

% Foulder for all source files recursively
addpath(genpath([file_path 'src']));

% Foulder for all external files recursively
addpath(genpath([file_path 'external']));

% Foulder for all source files recursively
%addpath(genpath([file_path 'test']));

end
