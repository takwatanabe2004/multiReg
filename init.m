% addpath(genpath('./convenience'))
% addpath(genpath('./util'))
% addpath(genpath('./gridsearch'))
rootdir=fileparts(mfilename('fullpath'));
addpath(genpath(rootdir))
rmpath([rootdir,'/.git'])


%-------------------------------------------------------------------------%
rmpath([rootdir,'/admm_scripts/admm_scripts_old'])

%-------------------------------------------------------------------------%
% - remove directory not of interest (02/27/2014)
%-------------------------------------------------------------------------%
%%
format compact