% addpath(genpath('./convenience'))
% addpath(genpath('./util'))
% addpath(genpath('./gridsearch'))
rootdir=fileparts(mfilename('fullpath'));
addpath(genpath(rootdir))
rmpath([rootdir,'/.git'])


%-------------------------------------------------------------------------%
rmpath([rootdir,'/admm_scripts/admm_scripts_binary_classif/old_connectomes'])
rmpath([rootdir,'/admm_scripts/old'])
rmpath([rootdir,'/simulation/mni_brain_slice/old_stuffs'])
rmpath([rootdir,'/data_local/_old_stuffs'])
rmpath([rootdir,'/proto/adjmat_kronless'])
rmpath([rootdir,'/proto/diffmat_kronless'])
rmpath(genpath([rootdir,'/proto/diffmat_newkron']))
rmpath(genpath([rootdir,'/data_local/_from_other_repos']))
%-------------------------------------------------------------------------%
% - remove directory not of interest (02/27/2014)
%-------------------------------------------------------------------------%
%%
format compact