% sim_plot_gridsearch_result.m (01/07/2014)
% - superseded by sim_plot_gridsearch_result_ident.m...
%%
error('superseded by sim_plot_gridsearch_result_ident.m')
clear
purge
cmap=[0.48,0.95];
gridname='grid2';
ntr=100;
% ntr=200;
snr=0.6;
%%
datapath=[gridname,'_snr',num2str(snr),'_enet_ntr',num2str(ntr),'.mat']
enet=load(datapath)

figure,imexpb
subplot(121),imagesc(log2(enet.gamgrid),log2(enet.lamgrid),enet.accuracy'),colorbar,impixelinfo
caxis(cmap)
subplot(122),imagesc(log2(enet.gamgrid),log2(enet.lamgrid),log10(enet.nnz_array')),
colorbar,
impixelinfo
drawnow
%%
datapath=[gridname,'_snr',num2str(snr),'_gnet_ntr',num2str(ntr),'.mat']
gnet=load(datapath)

figure,imexpb
subplot(121),imagesc(log2(gnet.gamgrid),log2(gnet.lamgrid),gnet.accuracy'),colorbar,impixelinfo
caxis(cmap)
subplot(122),imagesc(log2(gnet.gamgrid),log2(gnet.lamgrid),log10(gnet.nnz_array')),
colorbar,
impixelinfo
drawnow
%%
datapath=[gridname,'_snr',num2str(snr),'_flas_ntr',num2str(ntr),'.mat']
flas=load(datapath)
figure,imexpb
subplot(121),imagesc(log2(flas.gamgrid),log2(flas.lamgrid),flas.accuracy'),colorbar,impixelinfo
caxis(cmap)
subplot(122),imagesc(log2(flas.gamgrid),log2(flas.lamgrid),log10(flas.nnz_array')),
colorbar,
impixelinfo
drawnow
%% find peak accuracy
% enet
disp('---------------------')
[~,j]=max(enet.accuracy(:));
[igam,ilam]=ind2sub([length(enet.gamgrid),length(enet.lamgrid)],j);
best_enet.accuracy = enet.accuracy(igam,ilam);
best_enet.nnz      = enet.nnz_array(igam,ilam);
best_enet.lam = enet.lamgrid(ilam);
best_enet.gam = enet.gamgrid(igam);
best_enet.log2lam = log2(best_enet.lam);
best_enet.log2gam = log2(best_enet.gam)


% gnet
[~,j]=max(gnet.accuracy(:));
[igam,ilam]=ind2sub([length(gnet.gamgrid),length(gnet.lamgrid)],j);
best_gnet.accuracy = gnet.accuracy(igam,ilam);
best_gnet.nnz      = gnet.nnz_array(igam,ilam);
best_gnet.lam = gnet.lamgrid(ilam);
best_gnet.gam = gnet.gamgrid(igam);
best_gnet.log2lam = log2(best_gnet.lam);
best_gnet.log2gam = log2(best_gnet.gam)


% flas
[~,j]=max(flas.accuracy(:));
[igam,ilam]=ind2sub([length(flas.gamgrid),length(flas.lamgrid)],j);
best_flas.accuracy = flas.accuracy(igam,ilam);
best_flas.nnz      = flas.nnz_array(igam,ilam);
best_flas.lam = flas.lamgrid(ilam);
best_flas.gam = flas.gamgrid(igam);
best_flas.log2lam = log2(best_flas.lam);
best_flas.log2gam = log2(best_flas.gam)
%% scatter plot of accuracy vs sparsity
% YLIM=[0.5,1];
% % purge
% figure,imexpb
% subplot(131),scatter(enet.nnz_array(:),enet.accuracy(:)),ylim(YLIM)
% subplot(132),scatter(gnet.nnz_array(:),gnet.accuracy(:)),ylim(YLIM)
% subplot(133),scatter(flas.nnz_array(:),flas.accuracy(:)),ylim(YLIM)
% 
% figure,imexpb
% subplot(131),scatter(log10(enet.nnz_array(:)),enet.accuracy(:)),ylim(YLIM)
% subplot(132),scatter(log10(gnet.nnz_array(:)),gnet.accuracy(:)),ylim(YLIM)
% subplot(133),scatter(log10(flas.nnz_array(:)),flas.accuracy(:)),ylim(YLIM)