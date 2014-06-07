%% may30_spectr_rad_flasso2.m
% (05/30/2014)
%*************************************************************************%
% same as may30_spectr_rad_flasso.m, but add identity matrix on top of diffmat
%=========================================================================%
% seems like circulant fused lasso diffmat has spectral norm of sqrt(d^2),
% where d is the dimension
%-------------------------------------------------------------------------%
% eg: 1d -> spectral norm = sqrt(4)
%     2d -> spectral norm = sqrt(8) = 2sqrt(2)
%-------------------------------------------------------------------------%
% For 2d case, see 2014 Y. Ouyang - AN ACCELERATED LINEARIZED ADMM.pdf, 
% Sec3.2, where there's an exmaple of 2d TV...sqrt(8) is mentioned there, 
% with reference to 2004 A. Chambolle 
%
% Also see ~~2011 Chambolle, Pock - A first-order primal-dual algorithm for...
% pg 8, non journal version, for this statement
%=========================================================================%
%%
clear
purge
%%
disp('------- 1d -----------')
NSIZE1d=500;
C1d=tak_diffmat_1d(NSIZE1d,1);
I1d=speye( prod(NSIZE1d) );
tic
sv_I1d=svds(I1d,1)^2
sv_C1d=svds(C1d,1)^2
sv_F1d=svds([I1d; C1d],1)^2
toc
% return
%%
disp('------- 2d -----------')
NSIZE2d=[100,100];
I2d=speye( prod(NSIZE2d) );
C2d=tak_diffmat_2d(NSIZE2d,1);
tic
sv_I2d=svds(I2d,1)^2
sv_C2d=svds(C2d,1)^2
sv_F2d=svds([I2d; C2d],1)^2
toc
% return
%%
disp('------- 3d -----------')
NSIZE3d=[33,33,33];
I3d=speye( prod(NSIZE3d) );
C3d=tak_diffmat_3d(NSIZE3d,1);
tic
sv_I3d=svds(I3d,1)^2
sv_C3d=svds(C3d,1)^2
sv_F3d=svds([I3d; C3d],1)^2
toc
% return
%%
disp('------- 4d -----------')
NSIZE4d=[15,15,15,15];
I4d=speye( prod(NSIZE4d) );
C4d=tak_diffmat_4d(NSIZE4d,1);
tic
sv_I4d=svds(I4d,1)^2
sv_C4d=svds(C4d,1)^2
sv_F4d=svds([I4d; C4d],1)^2
toc
% return