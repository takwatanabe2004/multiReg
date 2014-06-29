%% june28_kronless_diffmat_6d
% (06/28/2014)
%=========================================================================%
% - Comments
%=========================================================================%
%%
clear all;
purge
%% 2d trial
NSIZE=repmat([7,5,6]*2,[1,2]);
% NSIZE=[4,3,2,4,3,2];
p=prod(NSIZE)

%=========================================================================%
% circulant case
%=========================================================================%
% % kron version
% tic
% [Ccirc,Cdim] = tak_diffmat_6d(NSIZE,1);
% Ccirc_x=Cdim.CX1;
% Ccirc_y=Cdim.CY1;
% Ccirc_z=Cdim.CZ1;
% Ccirc_t=Cdim.CX2;
% Ccirc_u=Cdim.CY2;
% Ccirc_v=Cdim.CZ2;
% toc
% % return
% % nonkron version
% tic
% [Ccirc2,Ccirc_x2,Ccirc_y2,Ccirc_z2,Ccirc_t2,Ccirc_u2,Ccirc_v2]=...
%     tak_diffmat_6d_nokron(NSIZE,1);
% toc
% 
% %-------------------------------------------------------------------------%
% % sanity checks via figures and "isequal"
% %-------------------------------------------------------------------------%
% if p<5000
% %     imcovvl(Ccirc_u)
% %     imcovvl(Ccirc_u2),caxis([-1,1])
% %     imcovvl(Ccirc_u-Ccirc_u2),caxis([-1,1])
% %     return
%     figure,imexp,
%     subplot(141),imcov(Ccirc_x),subplot(142),imcov(Ccirc_y),
%     subplot(143),imcov(Ccirc_z),subplot(144),imcov(Ccirc_t),
%     figure,imexp,
%     subplot(141),imcov(Ccirc_x2),subplot(142),imcov(Ccirc_y2),
%     subplot(143),imcov(Ccirc_z2),subplot(144),imcov(Ccirc_t2)
% else
% %     figure,imexp,subplot(121),tspy(Ccirc_x),subplot(122),tspy(Ccirc_y)
% %     figure,imexp,subplot(121),tspy(Ccirc_x2),subplot(122),tspy(Ccirc_y2)
% end
% if isequal(Ccirc,Ccirc2) && isequal(Ccirc_x,Ccirc_x2) && isequal(Ccirc_y,Ccirc_y2)
%     disp('success!!!')
% else
%     error('...')
% end
% Whos Ccirc
% % return

%=========================================================================%
% non-circulant case
%=========================================================================%
% kron version
tic
[C,Cdim] = tak_diffmat_6d(NSIZE,0);
Cx=Cdim.CX1;
Cy=Cdim.CY1;
Cz=Cdim.CZ1;
Ct=Cdim.CX2;
Cu=Cdim.CY2;
Cv=Cdim.CZ2;
toc
% return
% nonkron version
tic
[C2,Cx2,Cy2,Cz2,Ct2,Cu2,Cv2]=tak_diffmat_6d_nokron(NSIZE,0);
toc

%-------------------------------------------------------------------------%
% sanity checks via figures and "isequal"
%-------------------------------------------------------------------------%
if p<5000
%     imcovvl(Cv)
%     imcovvl(Cv2),caxis([-1,1])
%     imcovvl(Cv-Cv2),caxis([-1,1])
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

