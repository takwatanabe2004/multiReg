%% july05_test_1d_patch_MTL.m
% (07/05/2014)
%=========================================================================%
% - Test out tak_sim_groundTruthWeight_1d_patch_MTL.m
%=========================================================================%
%%
clear all;
purge

p = 1500;
q = 50; 

[W,idx_supp] = tak_sim_weight_1d_patch_MTL(p,q);
tplott(idx_supp)
figure,imexpl,imagesc(W),impixelinfo,drawnow

figure,imexpb
subplot(141),tplot(W(:,1))
subplot(142),tplot(W(:,2))
subplot(143),tplot(W(:,3))
subplot(144),tplot(W(:,4))

figure,imexpb
plot(W(:,1:8),'linewidth',3),grid on