%% june28_augMat_and_maskMat_looplessFunc
% (06/28/2014)
%=========================================================================%
% - attempt to replicate save_augMat_and_maskMat.m w/o looping
% - same as june28_augMat_and_maskMat_loopless.m, but using a function
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
[A2,b2]=tak_get_augMat_maskMat(ARRAYSIZE,coord.slex);
% Whos
% return
if isequal(A,A2)
    disp('******** SUCCESS!!!!!! ************')
else
    error('i suck...')
end

%% test if i can get the same FL&GN penalty values using the above A & b
Ccirc=tak_diffmat(ARRAYSIZE,1);
w = randn(p,1);

PEN.FL_brute = norm(C*w,1);
PEN.FL_circ1 = norm(b.*(Ccirc*(A*w)),1);
PEN.FL_circ2 = norm(b2.*(Ccirc*(A2*w)),1);

PEN.GN_brute = norm(C*w,2)^2;
PEN.GN_circ1 = norm(b.*(Ccirc*(A*w)),2)^2;
PEN.GN_circ2 = norm(b2.*(Ccirc*(A2*w)),2)^2;

PEN
if isequal(PEN.FL_circ1,PEN.FL_circ2) && isequal(PEN.GN_circ1,PEN.GN_circ2)
    disp('******** SUCCESS!!!!!! ************')
else
    error('i suck...')
end
