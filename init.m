% addpath(genpath('./convenience'))
% addpath(genpath('./util'))
% addpath(genpath('./gridsearch'))
rootdir=fileparts(mfilename('fullpath'));
addpath(genpath(rootdir))
rmpath([rootdir,'/.git'])


%-------------------------------------------------------------------------%
rmpath([rootdir,'/admm_scripts/admm_scripts_binary_classif/old_connectomes'])
rmpath([rootdir,'/simulation/15by15/old_stuffs'])
rmpath([rootdir,'/data_local/_old_stuffs'])
%-------------------------------------------------------------------------%
% - remove directory not of interest (02/27/2014)
%-------------------------------------------------------------------------%
%%
format compact