% Simulation_EstimParam_check (01/07/2014)
% - sanity check on reproducing Simulation_EstimParam.m
% - SUCCESS!
%%
clear
purge
load rcorr_design_censor X y
load yeo_info347 roiMNI
slab=load('C:\Users\takanori\Desktop\structSVM\simulation_new\Results.mat');
% Whos
z_used=slab.coord_used(1,3);
nodemask=(z_used==roiMNI(:,3));
coord_used=roiMNI(nodemask,:);

isequal(nodemask,slab.nodemask)
isequal(coord_used,slab.coord_used)
%%
nNodes=sum(nodemask);

edgemask = zeros(347);

for i=1:347
    for j=i:347
        if nodemask(i) & nodemask(j)
            edgemask(i,j)=1;
        end
    end
end
isequal(edgemask,slab.edgemask)

% symmetrize
edgemask = edgemask + edgemask';
% imedgel(slab.edgemask)
% imedgel(edgemask)
edgemaskflat=tak_dvec(edgemask);

idx_HC=(y==-1);
X_HC=X(idx_HC,:);

X_HC_mini = X_HC(:,logical(edgemaskflat));
X_HC_mini_fisher = tak_fisher_transformation(X_HC_mini);

ConnMeanVec=mean(X_HC_mini_fisher)';
ConnVarVec=var(X_HC_mini_fisher,1)'; % normalize by 1/N...MLE

ConnMean=tak_dvecinv(ConnMeanVec,0);
ConnVar=tak_dvecinv(ConnVarVec,0);

figure,imexp
subplot(121),imcov(ConnMean)
subplot(122),imcov(ConnVar)

figure,imexp
subplot(121),imcov(slab.ConnMean+slab.ConnMean')
subplot(122),imcov(slab.ConnVar+slab.ConnVar')

figure,imexp
subplot(121),imcov(abs(ConnMean-slab.ConnMean-slab.ConnMean'))
subplot(122),imcov(abs(ConnVar - slab.ConnVar-slab.ConnVar'))