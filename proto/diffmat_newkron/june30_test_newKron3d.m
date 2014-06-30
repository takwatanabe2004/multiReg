%% june30_test_newKron6d
% (06/30/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge
%% 3d trial
NSIZE=[250 250 25];
% NSIZE=repmat([7,8,11]*1,[1,2]);
% NSIZE=repmat([11,11,8]*1,[1,2]);
p=prod(NSIZE)

tic
[C,CX,CY,CZ]=tak_diffmat_3d(NSIZE);
toc

tic
[C2,CX2,CY2,CZ2]=tak_diffmat_3d_newkron(NSIZE);
toc

isequal(C,C2)
%%
w=randn(p,1);
isequal(CX*w,CX2*w)
isequal(CY*w,CY2*w)
isequal(CZ*w,CZ2*w)
% C2*w