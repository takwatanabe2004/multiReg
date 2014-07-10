%% july07_save_correlatedSmoothOutput1d.m.m
% (07/07/2014)
%=========================================================================%
% - Test script for my new function for adding patches
% (tak_sim_assignPatch1d.m, tak_sim_assignPulse1d.m)
% - The approach here is "dictionary method"...
%-------------------------------------------------------------------------%
% - the code looks unholy as of now...try making a funciton to modularize
%   this monstrosity
%=========================================================================%
%%
clear all;
purge

% randn('state',0)
% rand('state',0)

p = 1e3;
q = 15;

%=========================================================================%
% numClusters: determines the number of "clusters" in the set of
%              output weight vectors.  Each "clustesr" of weight vectors
%              will receive (weighted) basis elemetns from its own dictnoary
%   d: # atoms per dictinoary
%   nsubset: # dictinoary atoms a weight vector receives.  If this equals d,
%            each cluster should have the exact same shared supprot
%=========================================================================%
numClusters = 3; % <- determines # of dictionaries
d = 5; % <- # atoms per dictionary
nsubset = 5; % <- # dictionary atoms a weight vector receives
%% create (q x q) positive definite matrix
S = zeros(q,q);


idxCluster = cell(numClusters,1);
offset = 0;
for iClus = 1:numClusters
    %%%% for numClusters clusters of correlated features %%%%%%
    len=round(q/numClusters);
    
    if iClus==numClusters
        idxCluster{iClus} = offset+1:q;
    else
        idxCluster{iClus} = offset + (1:len);
        offset = offset + len;
    end
end
% return
%%%%%%%%%%%%%%%%%%%%%%%
%% collect "dictionary" atoms
Dlist = cell(numClusters,1);
for ind_c = 1:numClusters    
    D = zeros(p,d); % dictionary matrix
    for id=1:d
%         nPatches = 3;
        nPatches = randsample(3,1) + 2; % {3,...,5}
        pulseLenBase = round(p/500)*(randsample(2,1)+2);
        D(:,id) = tak_sim_assignPulse1d(p,nPatches,pulseLenBase)/50;
    end
    Dlist{ind_c}=D;
end
figure,imexp
for ii=1:min(numClusters,6)
    subplot(2,3,ii),imagesc(Dlist{ii}),impixelinfo,title('D = dictionary matrix')
end
% return
%% construct weight matrix cluster at a time
W = zeros(p,q);

% for icluster = 1:numClusters
%     for ii=1:length(idxCluster{icluster})
%         ind_set = idxCluster{icluster}(ii);
% %         randVec = 0*tak_sample_signed_unif([3,5],d) + 1*randn(d,1)
%         randVec = repmat(tak_sample_signed_unif([2,5],1),[d,1]) + 0*randn(d,1)
% %         randVec = sign(randn(d,1)).*repmat(tak_sample_signed_unif([3,5],1),[d,1]) + 1*randn(d,1)
%         tmp = Dlist{icluster}*randVec;
%         W(:,ind_set) = tmp;
%     end
% end

for icluster = 1:numClusters
    sizeCluster = length(idxCluster{icluster});
    %%
%     mixingMatrix = randn(d,sizeCluster)
    mixerBase = tak_sample_signed_unif([3,5],d)
%     mixingMatrix = repmat(mixerBase,[1,sizeCluster]) + 1*randn(d,sizeCluster)
    mixingMatrix = mixerBase*sign(randn(1,sizeCluster)) + 4*randn(d,sizeCluster)
    %%
%     W(:,idxCluster{icluster}) = Dlist{icluster}*mixingMatrix;

    % pick out subset of the mixer
    for iii=1:sizeCluster
        idx_subset = randsample(d,nsubset); % <- select nsubset dictionary atoms
        W(:,idxCluster{icluster}(iii)) = Dlist{icluster}(:,idx_subset)*mixingMatrix(idx_subset,iii); 
    end
end
% return
figure,imexpb,
subplot(121),imsupp(W)
subplot(122),imagesc(W),colormap jet,colorbar,impixelinfo
% return
figure,imexpb
subplot(121),imcov(corr(W)),title('corr(W)')
subplot(122),imcov(abs(corr(W))>0.3),title('corr(W)>0.3')

% figure,imexp
% for ii=1:min(q,12)
%     subplot(3,4,ii),tplot(W(:,ii))
% end
%% sample stuffs
n = 500;

%=========================================================================%
% iid normal features
%=========================================================================%
Xiid = randn(n,p);
sig=1;
noise=sig*randn(n,q);
Yiid = Xiid*W + noise;

%=========================================================================%
% what if features have structured distribution?
%=========================================================================%
Xar = tak_sample_AR1d(p,0.99,n);
Yar = Xar*W + noise;

figure,imexpb
subplot(121),imcov(corr(Yiid),1),  title('Correlated output')
subplot(122),imcov(corr(Yar),1),  title('Correlated output (AR features)')
figure,imexpb
subplot(121),imcov(abs(corr(Yiid))>0.3),  title('Correlated output >0.3')
subplot(122),imcov(abs(corr(Yar))>0.3),  title('Correlated output (AR features) >0.3')

tplottl(Yar)
%%
% purge
figure,imexp
subplot(311),tplot(W(:,idxCluster{1}))
subplot(312),tplot(W(:,idxCluster{2}))
subplot(313),tplot(W(:,idxCluster{3}))
%% save
rootdir=fileparts(mfilename('fullpath'));
mFileName=mfilename;
timeStamp=tak_timestamp;
% save([rootdir,'/correlatedWeight_1d_smooth2.mat'],'W','idxCluster','mFileName','timeStamp')
