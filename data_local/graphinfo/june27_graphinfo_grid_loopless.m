%% june27_graphinfo_grid_loopless.m
% (06/27/2014)
%=========================================================================%
% - attempt to replicate save_graphinfo_grid.m w/o looping
% - code mostly from june27_try_sim_347_graphinfo_2d_loopless.m, where I 
%   tested on a 4-d simulation connectome
%=========================================================================%
%%
clear
purge

rootdir = fileparts(mfilename('fullpath'));
parcellation = 'Grid1068'; % {'Grid326','Grid1068'}
%% 
dataPath=[rootdir, '/graph_info_',parcellation,'.mat'];
dataVars={'adjmat', 'C', 'coord', 'roiMNI', 'roiMNI_flipx','timeTotal'};
load(dataPath,dataVars{:})

timeTotal

dd=prod(coord.NSIZE); % #nodes in full FOV (includes "ghost" nodes)

%-------------------------------------------------------------------------%
% C6d=tak_diffmat_6d( [coord.NSIZE, coord.NSIZE],0);
% Lap6d=C6d'*C6d;
% adj6d = tak_inc2adj(C6d);
%-------------------------------------------------------------------------%
adj6d = tak_adjmat([coord.NSIZE, coord.NSIZE]);

% figure,imexpb
% subplot(131),tspy(C6d)
% subplot(132),tspy(Lap6d)
% subplot(133),tspy(adj6d)
% size(C6d,1)/size(C,1)

%-------------------------------------------------------------------------%
% create a dd^2-length mask-vector indicating 6-d coordinates included in
% the connectome
%-------------------------------------------------------------------------%
% mask_vec=false(dd^2,1);
% mask_vec(coord.slex)=true;
% Whos

%-------------------------------------------------------------------------%
% sample adjmat indices cooresponding to the above mask
%-------------------------------------------------------------------------%
adjmat_subsamp = adj6d(coord.slex,coord.slex);
C6d_subsamp = tak_adjmat2incmat(adjmat_subsamp);
Lap6d_subsamp=C6d_subsamp'*C6d_subsamp;
% figure,imexpb
% subplot(131),tspy(C6d_subsamp)
% subplot(132),tspy(Lap6d_subsamp)
% subplot(133),tspy(adjmat_subsamp)

%=========================================================================%
% compare with the brute-force method
%=========================================================================%
C6d_brute = tak_adjmat2incmat(adjmat);
Lap6d_brute=C6d_brute'*C6d_brute;
% figure,imexpb
% subplot(131),tspy(C6d_brute)
% subplot(132),tspy(Lap6d_brute)
% subplot(133),tspy(adjmat)

isequal(C6d_brute,C6d_subsamp)
isequal(adjmat,adjmat_subsamp)