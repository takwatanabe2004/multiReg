clear all
% purge

penalty='isoTV'; % {'enet','gnet','flas','isoTV'}

gridname = 'grid_may11';

dataPath=[get_rootdir,'/admm_scripts/gridsearch/',gridname,'_',...
    penalty,'.mat'];
%% load data
inputVars={'accuracy','gamgrid','lamgrid','accuracy','mean_nnz_array'};
gridresult=load(dataPath,inputVars{:});
% return
%%
% set figure options
fsize = 16; % fontsize 
fontOpt={'fontsize',fsize,'fontweight','b'};
xstr='log_2(\gamma)';
ystr='log_2(\lambda)';

axesOpt={'Fontweight','b','fontsize',fsize};

% colorbar range for performance
crangeAcc = [0.5,0.7]; % performance
crangeNNZ = [0 log10(60e3)]; % nnz
%%
gam=log2(gridresult.gamgrid);
lam=log2(gridresult.lamgrid);
figure,imexp
subplot(121),imagesc(gam,lam,gridresult.accuracy'),impixelinfo
xlabel(xstr,fontOpt{:}),ylabel(ystr,fontOpt{:}),colorbar
% caxis(crangeAcc)
set(gca,fontOpt{:})
title([penalty,'(accuracy)'])
subplot(122),imagesc(gam,lam,log10(gridresult.mean_nnz_array)'),impixelinfo
xlabel(xstr,fontOpt{:}),ylabel(ystr,fontOpt{:}),colorbar,caxis(crangeNNZ),
set(gca,fontOpt{:})
title([penalty,'(mean sparsity)'])
