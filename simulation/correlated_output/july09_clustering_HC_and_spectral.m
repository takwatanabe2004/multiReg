%% july09_clustering_HC_and_spectral.m
% (07/09/2014)
%=========================================================================%
% - Play with clustering algorithm to see if I can get the cluster groups
%   out of the output correlation matrix
%=========================================================================%
%%
clear all;
purge

randn('state',2)
rand('state',0)

rootdir = fileparts(mfilename('fullpath'));

load([rootdir,'/correlatedWeight_1d_patch3.mat'],'W','idxCluster')
% load([rootdir,'/correlatedWeight_1d_smooth2.mat'],'W','idxCluster')
%%
[p,q]=size(W);
%% make data
n=800;
ntest=100;

sig=5;
X=randn(n,p);
Xtest=randn(ntest,p);
% X = tak_sample_AR1d(p,0.8,n);
% Xtest = tak_sample_AR1d(p,0.8,ntest);

Y = X*W + sig*randn(n,q);

Ycorr = corr(Y);
%%
% play with spectral & hierarchical clustering
%=========================================================================%
% play with spectral clustering
%=========================================================================%
idx_perm=randsample(q,q);
Wtmp = Ycorr(idx_perm,idx_perm);
W= abs(Wtmp)>0.2;
% imcovvl(W)
% return
D = diag(sum(W,2));
L = D-W;
% L = eye(q)-D\W; % <- normalized graph laplacian
% [U,V]=eig(L);
[U,V]=eig(L);

K=3;
yy = U(:,1:K);

[idx,C]=kmeans(yy,K);
[~,idxsrt]=sort(idx)


figure,imexp
subplot(131), imcov(Wtmp)
subplot(132), imcov(Wtmp(idxsrt,idxsrt))
subplot(133), imcov(Ycorr)
% subplot(133), tstem(V)


%=========================================================================%
% play with (Agglomerative) hierarchical clustering
%=========================================================================%
Z=linkage(Wtmp.^2,'ward')
figure,imexpl,[~,~,idx3]=dendrogram(Z)

imcovvl(Wtmp)
imcovvl(Wtmp(idx3,idx3))


% imcovvl(Ycorr)