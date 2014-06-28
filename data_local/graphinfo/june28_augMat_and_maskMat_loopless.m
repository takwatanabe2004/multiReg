%% june28_augMat_and_maskMat_loopless
% (06/28/2014)
%=========================================================================%
% - attempt to replicate save_augMat_and_maskMat.m w/o looping
%=========================================================================%
%%
clear all;
purge

rootdir = fileparts(mfilename('fullpath'));
parcellation = 'Grid326'; % {'Grid326','Grid1068'}
%% load data from the save_augMat_and_maskMat.m function
dataPath = [rootdir,'/augmat_maskmat_',parcellation,'.mat'];
dataVars={'b', 'A', 'timeStamp',  'mFileName'};
load(dataPath,dataVars{:})
%% create diffmat looplessly (see june27_graphinfo_grid_looplessFunc.m)
% get coord info
load([rootdir, '/graph_info_',parcellation,'.mat'],'coord','C')

ARRAYSIZE=[coord.NSIZE,coord.NSIZE];

% adjmat = tak_adjmat_subsampled(ARRAYSIZE,coord.slex);
% C = tak_adjmat2incmat(adjmat);

dd=prod(coord.NSIZE); % #nodes including the "ghost" nodes
pp=dd^2; % #edges in ENTIRE 6-d space (including upper-triangular part of corrmat)

%-------------------------------------------------------------------------%
% p = #actual edges in the FC data 
% (consider this as the "subsampled" point from the full pp-cordinates)
%-------------------------------------------------------------------------%
[e,p]=size(C); 
%% try to replicate augmat looplessly
A2=sparse(coord.slex, 1:p, 1, pp,p);
Whos
return
if isequal(A,A2)
    disp('******** SUCCESS!!!!!! ************')
else
    error('i suck...')
end
% w=randn(p,1);
% isequal(A*w,A2*w)
%% try to replicate mask vector
disp('naw, just use the same method as before (was loopless anyways')
% NDIM = length(ARRAYSIZE);
% ee = NDIM * pp; 
% 
% C_circ = tak_diffmat(ARRAYSIZE,1);
% 
% b2 = false(ee,1);
% b3 = sparse(b);
% Whos b b2 b3