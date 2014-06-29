%% june28_kronless_adjmat_2d.m
% (06/28/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge
%% 2d trial
NSIZE=[2000,1500];
% NSIZE=[9,4];
p=prod(NSIZE);
%% non-circulant case
% % kron version
% [A,Ax,Ay] = tak_adjmat_2d(NSIZE,0);
% figure,imexp,subplot(131),tspy(A),subplot(132),tspy(Ax),subplot(133),tspy(Ay)
% % return
% 
% % nonkron version
% [A2,Ax2,Ay2] = tak_adjmat_2d_nokron(NSIZE,0);
% figure,imexp,subplot(131),tspy(A2),subplot(132),tspy(Ax2),subplot(133),tspy(Ay2)

%-------------------------------------------------------------------------%
% time battle
%-------------------------------------------------------------------------%
tic
A = tak_adjmat_2d(NSIZE,0);
t_Kron=toc

tic
A2 = tak_adjmat_2d_nokron(NSIZE,0);
t_nonKron=toc
if ~isequal(A,A2)
    error('debug!')
else
    disp('--- SUCCESS !!! ---')
end
return

%% circulant case
% kron version
tic
[Ccirc,Ccirc_x,Ccirc_y] = tak_diffmat_2d(NSIZE,1);
toc

% nonkron version
tic
[Ccirc2,Ccirc_x2,Ccirc_y2]=tak_diffmat_2d_nokron(NSIZE,1);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
    figure,imexp,subplot(121),imcov(Ccirc_x),subplot(122),imcov(Ccirc_y)
    figure,imexp,subplot(121),imcov(Ccirc_x2),subplot(122),imcov(Ccirc_y2)
else
%     figure,imexp,subplot(121),tspy(Ccirc_x),subplot(122),tspy(Ccirc_y)
%     figure,imexp,subplot(121),tspy(Ccirc_x2),subplot(122),tspy(Ccirc_y2)
end
if isequal(Ccirc,Ccirc2) && isequal(Ccirc_x,Ccirc_x2) && isequal(Ccirc_y,Ccirc_y2)
    disp('success!!!')
else
    error('...')
end
Whos C Ccirc