%% may30_spectr_rad_flasso.m
% (05/30/2014)
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
NSIZE1d=100;
C1d=full(tak_diffmat_1d(NSIZE1d,1));
L1d=C1d'*C1d;
tic
sv_1d=svds(C1d,1)^2
toc
% tic
% eigs(sparse(C1d))
% toc
[U1,V1]=eig(L1d);
U1,diag(V1)'
% %%
% disp('------- 2d -----------')
% NSIZE2d=[100,100];
% C2d=tak_diffmat_2d(NSIZE2d,1);
% tic
% sv_2d=svds(C2d,1)^2
% toc
% %%
% disp('------- 3d -----------')
% NSIZE3d=[33,33,33];
% C3d=tak_diffmat_3d(NSIZE3d,1);
% tic
% sv_3d=svds(C3d,1)^2
% toc
% %%
% disp('------- 4d -----------')
% NSIZE4d=[15,15,15,15];
% C4d=tak_diffmat_4d(NSIZE4d,1);
% tic
% sv_4d=svds(C4d,1)^2
% toc