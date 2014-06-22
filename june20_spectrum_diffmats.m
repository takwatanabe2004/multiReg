%=========================================================================%
% study the spectrum behavior of the difference matrices
%=========================================================================%
clear
purge

%% 1d
p=50;
C_1d=tak_diffmat_1d(p,0);
L_1d = C_1d'*C_1d;
svd_1d=svds(L_1d,p);

C_circ1d=tak_diffmat_1d(p,1);
L_circ1d = C_circ1d'*C_circ1d;
svd_circ1d=svd(full(L_circ1d));
figure,imexpl,tplot2(svd_1d,svd_circ1d),
legend('1d','1d-circ')
% imcovvl(L_circ1d)
% imcovvl(L_1d)
%% 2d
NSIZE=[41,41];
C_2d=tak_diffmat_2d(NSIZE,0);
L_2d = C_2d'*C_2d;
tic
svd_2d=svd(full(L_2d));
toc

C_circ2d=tak_diffmat_2d(NSIZE,1);
L_circ2d = C_circ2d'*C_circ2d;
svd_circ2d=svd(full(L_circ2d));
figure,imexpl,tplot2(svd_2d,svd_circ2d)
legend('2d','2d-circ')
% imcovvl(L_2d)
% imcovvl(L_circ2d)
%% 2d
NSIZE=[18,11,11];
C_3d=tak_diffmat_3d(NSIZE,0);
L_3d = C_3d'*C_3d;
tic
svd_3d=svd(full(L_3d));
toc

C_circ3d=tak_diffmat_3d(NSIZE,1);
L_circ3d = C_circ3d'*C_circ3d;
svd_circ3d=svd(full(C_circ3d));
figure,imexpl,tplot2(svd_3d,svd_circ3d)
legend('3d','3d-circ')