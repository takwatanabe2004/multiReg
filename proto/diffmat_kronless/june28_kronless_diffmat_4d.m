%% june28_kronless_diffmat_4d
% (06/28/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge
%% 2d trial
NSIZE=[7,5,6,4]*8;
% NSIZE=[7,5,6,4];
p=prod(NSIZE)

%=========================================================================%
% circulant case
%=========================================================================%
% kron version
tic
[Ccirc,Ccirc_x,Ccirc_y,Ccirc_z,Ccirc_t] = tak_diffmat_4d(NSIZE,1);
toc
if p<5000
    figure,imexp,
    subplot(141),imcov(Ccirc_x),subplot(142),imcov(Ccirc_y),
    subplot(143),imcov(Ccirc_z),subplot(144),imcov(Ccirc_t),
end
% return
% nonkron version
tic
[Ccirc2,Ccirc_x2,Ccirc_y2,Ccirc_z2,Ccirc_t2]=tak_diffmat_4d_nokron(NSIZE,1);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
%     imcovvl(Ccirc_y)
%     imcovvl(Ccirc_y2),caxis([-1,1])
%     imcovvl(Ccirc_y-Ccirc_y2),caxis([-1,1])
%     return
    figure,imexp,
    subplot(141),imcov(Ccirc_x),subplot(142),imcov(Ccirc_y),
    subplot(143),imcov(Ccirc_z),subplot(144),imcov(Ccirc_t),
    figure,imexp,
    subplot(141),imcov(Ccirc_x2),subplot(142),imcov(Ccirc_y2),
    subplot(143),imcov(Ccirc_z2),subplot(144),imcov(Ccirc_t2)
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
[C,Cx,Cy,Cz,Ct] = tak_diffmat_4d(NSIZE,0);
toc
% return
% nonkron version
tic
[C2,Cx2,Cy2,Cz2,Ct2]=tak_diffmat_4d_nokron(NSIZE,0);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
%     imcovvl(Cx)
%     imcovvl(Cx2),caxis([-1,1])
%     imcovvl(Cx-Cx2),caxis([-1,1])
%     return
    figure,imexp,
    subplot(141),imcov(Cx),subplot(142),imcov(Cy),
    subplot(143),imcov(Cz),subplot(144),imcov(Ct),
    figure,imexp,
    subplot(141),imcov(Cx2),subplot(142),imcov(Cy2),
    subplot(143),imcov(Cz2),subplot(144),imcov(Ct2)
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

