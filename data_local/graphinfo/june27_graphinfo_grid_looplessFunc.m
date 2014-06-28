%% june27_graphinfo_grid_looplessFunc
% (06/27/2014)
%=========================================================================%
% - same as june27_graphinfo_grid_loopless.m, but using a function
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

%=========================================================================%
% loopless construction of the differencing matrix
%=========================================================================%
adjmat_subsamp = tak_adjmat_subsampled([coord.NSIZE, coord.NSIZE],coord.slex);
C6d_subsamp = tak_adjmat2incmat(adjmat_subsamp);
Lap6d_subsamp=C6d_subsamp'*C6d_subsamp;

%=========================================================================%
% compare with the brute-force loop method
%=========================================================================%
C6d_brute = tak_adjmat2incmat(adjmat);
Lap6d_brute=C6d_brute'*C6d_brute;
% figure,imexpb
% subplot(131),tspy(C6d_brute)
% subplot(132),tspy(Lap6d_brute)
% subplot(133),tspy(adjmat)

isequal(C6d_brute,C6d_subsamp)
isequal(adjmat,adjmat_subsamp)
isequal(Lap6d_brute,Lap6d_subsamp)