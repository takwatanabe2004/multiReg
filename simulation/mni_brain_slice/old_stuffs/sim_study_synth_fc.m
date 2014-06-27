% sim_study_synth_fc.m (01/07/2014)
% - create synthetic functional connectome for simulation experiment
%==========================================================================
%--------------------------------------------------------------------------
%%
clear all
purge
load('Results_tak.mat', 'ConnMean', 'ConnVar')
load('graph_info347_2d.mat', 'coord')
load('sim_anom_node_info.mat', 'anom_nodes')

% set SNR devel (same as cohen's d) http://en.wikipedia.org/wiki/Cohen%27s_d#Cohen.27s_d 
% http://www.3rs-reduction.co.uk/html/6__power_and_sample_size.html
snr=1;

% nHC = number of control samples
% nDS = number of patient samples
nHC = 500;
nDS=nHC;
randn('seed',0)
%%
ConnMeanVec = tak_dvec(ConnMean);
ConnVarVec  = tak_dvec(ConnVar);
p=length(ConnMeanVec);
% tplott(ConnMeanVec)
% tplott(ConnVarVec)
%% create -1 (HC) control subjects
X_HC=repmat(ConnMeanVec(:)',[nHC,1])+ ...
     repmat(sqrt(ConnVarVec(:))',[nHC,1]).*randn(nHC,p);
 % figure,imagesc(X_HC),drawnow,colorbar,axis on,impixelinfo

%--------------------------------------------------------------------------
% sanity check: 
% - if nHC is large, the sample mean and the true mean should come close
%--------------------------------------------------------------------------
% figure,imexpb
% tmp=tak_dvecinv(mean(X_HC),0);
% subplot(131),imcov(tmp)
% subplot(132),imcov(ConnMean)
% subplot(133),imcov(abs(tmp-ConnMean))
% norm(tmp-ConnMean)
% 
% tmp=tak_dvecinv(var(X_HC,1),0);
% figure,imexpb
% subplot(131),imcov(tmp)
% subplot(132),imcov(ConnVar)
% subplot(133),imcov(abs(tmp-ConnVar)),drawnow
% norm(tmp-ConnVar)

%==========================================================================
% inverse fisher tx to map to (-1,+1) correlation space
%==========================================================================
% figure,hist(X_HC(:),1e3)
X_HC = tak_fisher_itransformation(X_HC);
% figure,hist(X_HC(:),1e3)
%% create +1 (DS) disease subjects
X_DS=repmat(ConnMeanVec(:)',[nDS,1])+ ...
     repmat(sqrt(ConnVarVec(:))',[nDS,1]).*randn(nDS,p);

% before signal is added
X_DS_before=tak_fisher_itransformation(X_DS);

%==========================================================================
% add SIGNAL on connections between the two anomalous node clusters
%==========================================================================
for i=1:length(anom_nodes.idx_conn)
    % connection index
    idx=anom_nodes.idx_conn(i);
    
    SIGNAL=snr*sqrt(ConnVarVec(idx));
    
    X_DS(:,idx)=X_DS(:,idx) + repmat(SIGNAL,[nDS,1]);
end
X_DS = tak_fisher_itransformation(X_DS);
[min(X_DS(:)),max( X_DS(:))]
% figure,imexpb
% subplot(131),imagesc(X_DS),drawnow,colorbar,axis on,impixelinfo
% subplot(132),imagesc(X_DS_before),drawnow,colorbar,axis on,impixelinfo
% subplot(133),imagesc(X_DS-X_DS_before),drawnow,colorbar,axis on,impixelinfo
figure,imexpb,imagesc(X_DS-X_DS_before),drawnow,colorbar,axis on,impixelinfo
max(max(abs(X_DS-X_DS_before)))

figure,imexpb
subplot(131),tplot(X_DS(1,:))
subplot(132),tplot(X_DS_before(1,:))
subplot(133),tplot(X_DS(1,:)-X_DS_before(1,:))
drawnow
max( X_DS(:))
max( X_DS_before(:))
%% plot mean between HC and DS
meanHC=mean(X_HC);
meanDS=mean(X_DS);
meanHC_mat=tak_dvecinv(meanHC,0);
meanDS_mat=tak_dvecinv(meanDS,0);

figure,imexpb
subplot(131),tplot(meanDS)
subplot(132),tplot(meanHC)
subplot(133),tplot(meanDS-meanHC)

figure,imexpb
subplot(131),imcov(meanDS_mat)
subplot(132),imcov(meanHC_mat)
subplot(133),imcov(meanDS_mat-meanHC_mat)
 
figure,imexpb
subplot(131),imcov(ConnMean)

drawnow