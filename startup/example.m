%% Example startup file to get the MATLAB environment ready
%  this stuff probably should have been included as submodules... you live and learn!

% https://github.com/opencobra/cobratoolbox
addpath(genpath('~/dev/cobratoolbox'));

% CPLEX
addpath(genpath('/Applications/CPLEX_Studio129/cplex/matlab'));
changeCobraSolver('cplex_direct', 'LP');

% this repo -- for funtions in subfolders
addpath(genpath('~/dev/epfl-systems-biology'));

% https://github.com/EPFL-LCSB/matTFA
addpath(genpath('~/dev/matTFA'));
