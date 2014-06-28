%% save_multimodal_crosscorrmat
% (05/11/2014)
%=========================================================================%
%=========================================================================%
%%
clear
purge

X1=struct2array(load([get_rootdir,'/data_local/sMRI_design_censor.mat'],'X'));
X2=struct2array(load([get_rootdir,'/data_local/rcorr_design_censor.mat'],'X'));
[n,p1]=size(X1);
[~,p2]=size(X2);

fsave=true;

threshold = 0.3;

outputPath=[get_rootdir,'/data_local/multimodal_crosscorrmat',...
    num2str(threshold),'.mat']
% return
%%
% X1z=zscore(X1);
% meanX1z=mean(X1z);
% varX1z=std(X1z);
% 
% X2z=zscore(X2);
% meanX2z=mean(X2z);
% varX2z=std(X2z);

%-------------------------------------------------------------------------%
% for the sMRI data, voxels that are 0-valued for all 121 subjects have 
% 0 variance from the zscore standardization
%-------------------------------------------------------------------------%
% sum(sum(bsxfun(@eq,X1,0))==121)
% sum(varX1z==0)
%%
cross_corrmat = sparse(p1,p2);
tic
for i=1:p1
    if mod(i,500)==0, toc, disp([num2str(i),' out of ',num2str(p1)]), tic,end;
    
    cross_corr=corr(X1(:,i),X2);
    
    mask=abs(cross_corr)>threshold;
    cross_corrmat(i,mask) = cross_corr(mask);
    
%     cross_corr(mask)
end
% return

timeStamp=tak_timestamp;
mFileName=mfilename;
if fsave
    save(outputPath, 'cross_corrmat', 'timeStamp','mFileName')
end