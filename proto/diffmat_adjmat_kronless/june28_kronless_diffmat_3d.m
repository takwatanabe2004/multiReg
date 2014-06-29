%% june28_kronless_diffmat_3d
% (06/28/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge
%% 2d trial
NSIZE=[100,200,150];
% NSIZE=[8,3,4];
p=prod(NSIZE)

%=========================================================================%
% circulant case
%=========================================================================%
% kron version
tic
[Ccirc,Ccirc_x,Ccirc_y,Ccirc_z] = tak_diffmat_3d(NSIZE,1);
toc
if p<5000
    figure,imexp,
    subplot(131),imcov(Ccirc_x),
    subplot(132),imcov(Ccirc_y),
    subplot(133),imcov(Ccirc_z)
end
% return
% nonkron version
tic
[Ccirc2,Ccirc_x2,Ccirc_y2,Ccirc_z2]=tak_diffmat_3d_nokron(NSIZE,1);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
%     imcovvl(Ccirc_x)
%     imcovvl(Ccirc_x2),caxis([-1,1])
%     imcovvl(Ccirc_x-Ccirc_x2),caxis([-1,1])
%     return
    figure,imexp,
    subplot(131),imcov(Ccirc_x),subplot(132),imcov(Ccirc_y),
    subplot(133),imcov(Ccirc_z)
    figure,imexp,
    subplot(131),imcov(Ccirc_x2),subplot(132),imcov(Ccirc_y2),
    subplot(133),imcov(Ccirc_z2)
else
%     figure,imexp,subplot(121),tspy(Ccirc_x),subplot(122),tspy(Ccirc_y)
%     figure,imexp,subplot(121),tspy(Ccirc_x2),subplot(122),tspy(Ccirc_y2)
end
if isequal(Ccirc,Ccirc2) && isequal(Ccirc_x,Ccirc_x2) && isequal(Ccirc_y,Ccirc_y2)
    disp('success!!!')
else
    error('...')
end
Whos Ccirc
% return

%=========================================================================%
% non-circulant case
%=========================================================================%
% kron version
tic
[C,Cx,Cy,Cz] = tak_diffmat_3d(NSIZE,0);
toc
if p<5000
    figure,imexp,subplot(131),imcov(Cx),subplot(132),imcov(Cy)
                 subplot(133),imcov(Cz)
end
% return
% nonkron version
tic
[C2,Cx2,Cy2,Cz2]=tak_diffmat_3d_nokron(NSIZE,0);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
%     imcovvl(Cx)
%     imcovvl(Cx2),caxis([-1,1])
%     imcovvl(Cx-Cx2),caxis([-1,1])
    figure,imexp,subplot(131),imcov(Cx),subplot(132),imcov(Cy)
                 subplot(133),imcov(Cz)
    figure,imexp,subplot(131),imcov(Cx2),subplot(132),imcov(Cy2)
                 subplot(133),imcov(Cz2)
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

