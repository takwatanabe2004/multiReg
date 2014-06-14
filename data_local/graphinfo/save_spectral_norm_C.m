%% save_spectral_norm_C
% (06/10/2014)
%=========================================================================%
% - Precompute and save the spectral norm of the differencing matrix C
%  (required for setting the stepsizes of the C&P's prima-dual algorithm)
%=========================================================================%
%%
clear
purge

grid='Grid326'; % {'Grid326','Grid1068','WashU'}

outPath = [fileparts(mfilename('fullpath')),'/spectral_norm_C_',grid,'.mat']
outVars = {'spec_norm_C', 'mFileName', 'timeStamp'};
% return
%%
%=========================================================================%
% load finite differencing matrix
% (faster to compute the largest eigenvalue of the laplacian matrix here)
%=========================================================================%
load([get_rootdir,'/data_local/graphinfo/graph_info_',grid,'.mat'],'adjmat')
C = tak_adjmat2incmat(adjmat);

LAP = C'*C;
tic
% spec_norm_C = svds(C,1)^2 <- this approach is slower
spec_norm_C = sqrt(eigs(LAP,1));
toc
%%
mFileName=mfilename('fullpath');
timeStamp=tak_timestamp;
save(outPath,outVars{:})
%%

% %% compute spectral norm for step size
% %-------------------------------------------------------------------------%
% % spectral norm for step size of cppd alg (sigma*tau L^2 < 1 must be satisfied)
% %-------------------------------------------------------------------------%
% lam=11*rand
% gam=5*rand
% F = [lam*speye(p);gam*C];
% 
% 
% rad_C = svds(C,1)^2
% rad_gamC = svds(gam*C,1)^2
% 
% rad_F = svds(F,1)^2
% 
% % sqrt(rad_gamC/rad_C)
% %-------------------------------------------------------------------------%
% % LAP=C'*C;
% % degree=sum(abs(LAP)-diag(diag(LAP)));
% % tplott(sort(degree))
% % tplott(diag(LAP))
% %-------------------------------------------------------------------------%
% % svds(F,1)^2-svds(C,1)^2
% %%
% % s=svds(C,1)
% % svds(gam*C,1)
% % 
% % s*gam
% % (gam^2 * rad_C) + (lam^2)
% % 
% % sqrt(rad_F)
% % gam*svds(C,1) + lam