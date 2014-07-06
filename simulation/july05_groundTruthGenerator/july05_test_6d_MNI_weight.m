%% july05_test_6d_MNI_weight.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_weight_6d_MNI.m
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

[w,idx_supp,C] = tak_sim_weight_6d_MNI(clusterPath);

figure,imexp
subplot(121),tplot(idx_supp)
subplot(122),imconnEdge(idx_supp)
drawnow

figure,imexp
subplot(121),imconn(w,1)
subplot(122),imconnYeo(w,GRID,1,1)
drawnow