%% july09_save_correlatedPatchOutput2d.m
% (07/09/2014)
%=========================================================================%
% - Script built on top of july09_save_correlatedPatchOutput1d.m
%=========================================================================%
%%
clear all;
purge

% randn('state',0)
% rand('state',0)

nx=50;
ny=40;
NSIZE=[nx,ny];
p = prod(NSIZE);
q = 20;

%=========================================================================%
% nClusters: determines the number of "clusters" in the set of
%              output weight vectors.  Each "clustesr" of weight vectors
%              will receive (weighted) basis elemetns from its own dictnoary
%   natoms: # atoms per dictinoary
%   nsubset: # dictinoary atoms a weight vector receives.  If this equals d,
%            each cluster should have the exact same shared supprot
%=========================================================================%
% number of output block structure/clusters
nClusters = 4;

dict.natoms = 12; % <- # atoms per dictionary
dict.nsubset = 9; % <- # atoms per mixed signal (D *alpha, a
%% get indices of each cluster group (these are "output" indices)
idxCluster = cell(nClusters,1);
offset = 0;
for iClus = 1:nClusters
    %%%% for numClusters clusters of correlated features %%%%%%
    len=round(q/nClusters);
    
    if iClus==nClusters
        idxCluster{iClus} = offset+1:q;
    else
        idxCluster{iClus} = offset + (1:len);
        offset = offset + len;
    end
end
%% construct list of "dictionary" for each cluster group
dict.Dlist = cell(nClusters,1);
for ind_c = 1:nClusters    
    D = zeros(p,dict.natoms); % dictionary matrix
    %=====================================================================%
    % construct sparse dictionary matrix for each "cluster-group"
    % - dictionary atoms consist of smooth "pulse" signals
    %=====================================================================%
    for id=1:dict.natoms
        nPatches = 1;
%         nPatches = randsample(3,1) + 2; % {3,...,5}
        pulseLenBaseX = round(nx/25);
        pulseLenBaseY = round(ny/25);
        pulseLenBase = [pulseLenBaseX,pulseLenBaseY];
%         D(:,id) = tak_sim_assignPulse1d(p,nPatches,pulseLenBase)/50;
%         D(:,id) = tak_sim_assignPatch1d(p,nPatches,pulseLenBase)/50;
        D(:,id) = tak_sim_assignPatch2d(NSIZE,nPatches,pulseLenBase)/20;        
    end
%     keyboard
%     imcovvl(reshape(D(:,3),[NSIZE(1),NSIZE(2)]))
    dict.Dlist{ind_c}=D;
end
% return
figure,imexp
for ii=1:min(nClusters,6)
    subplot(2,3,ii),imagesc(dict.Dlist{ii}),impixelinfo,title('D = dictionary matrix')
end
% return
%% construct weight matrix cluster at a time
W = zeros(p,q);

%-------------------------------------------------------------------------%
% mixerCorr: dictates how "correlated" output structures are in each 
%            cluster blocks are
%-------------------------------------------------------------------------%
mixerCorr = 0.9; 
for icluster = 1:nClusters
    sizeCluster = length(idxCluster{icluster});
%     mixerCorr=-1*mixerCorr;
%     mixingMatrix = tak_sample_AR1d(sizeCluster, mixerCorr, dict.natoms);
    mixingMatrix = ones(dict.natoms,sizeCluster) .* ...
        repmat(sign(randn(1,sizeCluster)),[dict.natoms,1])...
        + .05*randn(dict.natoms,sizeCluster);
    
%     W(:,idxCluster{icluster}) = dict.Dlist{icluster}*mixingMatrix;
    % pick out subset of the mixer
    for iii=1:sizeCluster
        idx_subset = randsample(dict.natoms,dict.nsubset); % <- select nsubset dictionary atoms
        W(:,idxCluster{icluster}(iii)) = dict.Dlist{icluster}(:,idx_subset)*mixingMatrix(idx_subset,iii); 
    end
end

figure,imexpb
subplot(121),imagesc(W),colorbar,impixelinfo
subplot(122),imcov(corr(W),1),title('corr(W)'),colorbar off,colorbar Eastoutside


figure,imexpb,
subplot(121),imsupp(W)
subplot(122),imcov(abs(corr(W))>0.3),title('corr(W)>0.3'),colorbar off
%%
figure,imexp
subplot(141),imcov(reshape(W(:,idxCluster{1}(1)), [nx,ny]) ,1)
subplot(142),imcov(reshape(W(:,idxCluster{2}(1)), [nx,ny]) ,1)
subplot(143),imcov(reshape(W(:,idxCluster{3}(1)), [nx,ny]) ,1)
subplot(144),imcov(reshape(W(:,idxCluster{4}(1)), [nx,ny]) ,1)

figure,imexp
subplot(141),imcov(reshape(W(:,idxCluster{1}(2)), [nx,ny]) ,1)
subplot(142),imcov(reshape(W(:,idxCluster{2}(2)), [nx,ny]) ,1)
subplot(143),imcov(reshape(W(:,idxCluster{3}(2)), [nx,ny]) ,1)
subplot(144),imcov(reshape(W(:,idxCluster{4}(2)), [nx,ny]) ,1)
% return
%% sample stuffs
n = 500;

%=========================================================================%
% iid normal features
%=========================================================================%
Xiid = randn(n,p);
sig=0;
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

% tplottl(Yar)
%% save
rootdir=fileparts(mfilename('fullpath'));
mFileName=mfilename;
timeStamp=tak_timestamp;
% save([rootdir,'/correlatedWeight_july09_2d_patch2.mat'],'W','idxCluster','dict','nx','ny','mFileName','timeStamp')
