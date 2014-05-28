%% save_ttest_rank_rcorr
% (05/12/2014)
%=========================================================================%
% - save ttest ranking for the 10fold CV indices (for FC dataset)
%=========================================================================%
%%
clear
purge

fsave=true;
%% load data
load rcorr_design_censor.mat X y

[n,p]=size(X);

% scale to 0 mean, unit variance?
flag_scale = true;
if flag_scale
    X=zscore(X);
end

%==========================================================================
% load info for the indexing for the 10-fold-CV partitioning
%==========================================================================
load Overall Overall
%% output save path
if flag_scale
    outPath=[get_rootdir,'/data_local/rcorr_info/ttest_rank_idx_CV_zscaled.mat'];
else
    outPath=[get_rootdir,'/data_local/rcorr_info/ttest_rank_idx_CV.mat'];
end
outVars={'ttest_rank_idx_CV','mFileName','timeStamp'};
%%
for idxCV = 1:10
    idxCV
    %======================================================================
    % 10-fold-CV data partition
    %======================================================================
    mask_ts = Overall.CrossValidFold(:,idxCV);
    mask_tr = ~mask_ts;
    
    Xts = X(mask_ts,:);
    Xtr = X(mask_tr,:);

    yts = y(mask_ts);
    ytr = y(mask_tr);
    %======================================================================
    % get ttest rank
    %======================================================================
    ttest_rank_idx_CV{idxCV} = tak_ttest_prune(Xtr,ytr);
end
%% save
mFileName=mfilename;
timeStamp=tak_timestamp;

if fsave
    save(outPath,outVars{:})
end