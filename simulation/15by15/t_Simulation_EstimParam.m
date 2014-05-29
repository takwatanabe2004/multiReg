% t_Simulation_EstimParam.m (01/07/2014)
% - my version of creating the mean matrix and all that junk
% - sanity check conducted on Simulation_EstimParam_check.m
%==========================================================================
%--------------------------------------------------------------------------
%%
clear
purge

load rcorr_design_censor X y
load yeo_info347 roiMNI

% z-slice coordinate we use for the simulation
z_used = 18;
nodemask=(z_used==roiMNI(:,3));
coord_used=roiMNI(nodemask,:);

nNodes=sum(nodemask);
edgemask = false(347);

for i=1:347
    for j=i:347
        if nodemask(i) && nodemask(j)
            edgemask(i,j)=true;
            edgemask(j,i)=true;
        end
    end
end
% imedge(edgemask)
edgemaskVec=tak_dvec(edgemask);

%==========================================================================
% extract healthy subjects from dataset
%==========================================================================
idx_HC=(y==-1);
X_HC=X(idx_HC,:);

%==========================================================================
% extract edges corresponding to the z-slice of interest
%==========================================================================
X_HC_mini = X_HC(:,edgemaskVec);

%==========================================================================
% map correlation to (-inf,+inf) via fisher tx
%==========================================================================
X_HC_mini_fisher = tak_fisher_transformation(X_HC_mini);

ConnMeanVec=mean(X_HC_mini_fisher)';
ConnVarVec=var(X_HC_mini_fisher,1)'; % normalize by 1/N...MLE

ConnMean=tak_dvecinv(ConnMeanVec,0);
ConnVar=tak_dvecinv(ConnVarVec,0);

figure,imexp
subplot(121),imcov(ConnMean)
subplot(122),imcov(ConnVar)

%==========================================================================
% sanity check with the one daniel created
%==========================================================================
slab=load('C:\Users\takanori\Desktop\structSVM\simulation_new\Results.mat');
slab.ConnMean=slab.ConnMean+slab.ConnMean'; % symmetrize
slab.ConnVar=slab.ConnVar+slab.ConnVar'; % symmetrize

figure,imexp
subplot(121),imcov(slab.ConnMean)
subplot(122),imcov(slab.ConnVar)

figure,imexp
subplot(121),imcov(abs(ConnMean - slab.ConnMean))
subplot(122),imcov(abs(ConnVar - slab.ConnVar))

norm(ConnMean - slab.ConnMean)
norm(ConnVar - slab.ConnVar)
mFileName=mfilename;
timeStamp=tak_timestamp;
% save('Results_tak.mat','ConnMean','ConnVar','edgemask','edgemaskVec','coord_used','roiMNI','nodemask','mFileName','timeStamp')
