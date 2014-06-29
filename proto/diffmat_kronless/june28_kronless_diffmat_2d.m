%% june28_kronless_diffmat_2d
% (06/28/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge
%% 2d trial
% NSIZE=[2000,1500];
NSIZE=[20,15];
p=prod(NSIZE);

%=========================================================================%
% non-circulant case
%=========================================================================%
% kron version
tic
[C,Cx,Cy] = tak_diffmat_2d(NSIZE,0);
toc
% figure,imexp,subplot(121),imcov(Cx),subplot(122),imcov(Cy)

% nonkron version
tic
[C2,Cx2,Cy2]=tak_diffmat_2d_nokron(NSIZE,0);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
%     imcovvl(Cy)
%     imcovvl(Cy2),caxis([-1,1])
%     imcovvl(Cy-Cy2),caxis([-1,1])
    figure,imexp,subplot(121),imcov(Cx),subplot(122),imcov(Cy)
    figure,imexp,subplot(121),imcov(Cx2),subplot(122),imcov(Cy2)
else
%     figure,imexp,subplot(121),tspy(Ccirc_x),subplot(122),tspy(Ccirc_y)
%     figure,imexp,subplot(121),tspy(Ccirc_x2),subplot(122),tspy(Ccirc_y2)
end
if isequal(C,C2) && isequal(Cx,Cx2) && isequal(Cy,Cy2)
    disp('success!!!')
else
    error('...')
end
Whos C Ccirc

%=========================================================================%
% circulant case
%=========================================================================%
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