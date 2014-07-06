%% july05_test_6d_MNI_weight_MTL.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_6d_MNI_MTL.m
%=========================================================================%
%%
clear all;
purge

% rootdir = fileparts(mfilename('fullpath'));
%% load data
GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

%=========================================================================%
% anomalous edge index info
%=========================================================================% 
clusterName = 'anomCluster1';
clusterPath = [get_rootdir,'/simulation/july05_groundTruthGenerator/',clusterName,'_',GRID];

% # outputs
q = 20;
[W,idx_supp,C] = tak_sim_weight_6d_MNI_MTL(clusterPath,q);

figure,imexp
subplot(121),tplot(idx_supp)
subplot(122),imconnEdge(idx_supp)
drawnow

figure,imexpl,imagesc(W),impixelinfo,drawnow

figure,imexp
subplot(241),imconn(W(:,1),1)
subplot(242),imconn(W(:,2),1)
subplot(243),imconn(W(:,3),1)
subplot(244),imconn(W(:,4),1)
subplot(245),imconn(W(:,5),1)
subplot(246),imconn(W(:,6),1)
subplot(247),imconn(W(:,7),1)
subplot(248),imconn(W(:,8),1)

figure,imexp
subplot(241),imconnYeo(W(:,1),[],1)
subplot(242),imconnYeo(W(:,2),[],1)
subplot(243),imconnYeo(W(:,3),[],1)
subplot(244),imconnYeo(W(:,4),[],1)
subplot(245),imconnYeo(W(:,5),[],1)
subplot(246),imconnYeo(W(:,6),[],1)
subplot(247),imconnYeo(W(:,7),[],1)
subplot(248),imconnYeo(W(:,8),[],1)