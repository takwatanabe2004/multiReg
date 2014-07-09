%% july07_corrWeight_covMethod_fcn_patch1d.m.m
% (07/07/2014)
%=========================================================================%
% - Test script for my new function for adding patches
% (tak_sim_assignPatch1d.m, tak_sim_assignPulse1d.m)
% - The approach here is "covariance method"...
%=========================================================================%
%%
clear all;
purge

% randn('state',0)
% rand('state',0)

p = 5e3;
q = 20;
%% create (q x q) positive definite matrix
S = zeros(q,q);
%%%% for 3 clusters of correlated features %%%%%%
len1=round(q/3);
len2=round(q/3);
idx1= 1:len1;
idx2= 1+len1:len1+len2;
idx3= 1+len1+len2:q;
S(idx1,idx1)=0.9;
S(idx2,idx2)=0.8;
S(idx3,idx3)=-0.35; % <- when this is negative...S may not be PD
S = S - diag(diag(S)) + eye(q); % <- replace diagonal entries with 1
S=S+eye(q); % make S positive-definite
% [~,spectrum]=tak_eig(S); tstem(spectrum) % <- check PD condition of S
L=chol(S,'lower');
figure,imexpb
subplot(121),imcov(S)
subplot(122),imcov(L)
% return
%%%%%%%%%%%%%%%%%%%%%%%
%% 
W = zeros(p,q);
for iq=1:q
    nPatches = 2%randsample(4,1) + 4; % {5,...,7}
%     W(:,iq) = tak_sim_assignPatch1d(p,nPatches,round(p/200));
    W(:,iq) = tak_sim_assignPulse1d(p,nPatches,round(p/200));
end
Wnew = W*L;
figure,imexp
for ii=1:min(q,10)
    subplot(3,4,ii),tplot(W(:,ii))
end
figure,imexp
for ii=1:min(q,10)
    subplot(3,4,ii),tplot(Wnew(:,ii))
end
% return
figure,imexp
subplot(221),imcov(corr(W),1), title('corr(W)')
subplot(223),imagesc(W),caxis(max(abs(caxis))*[-1,1]); title('W')
subplot(222),imcov(corr(Wnew),1), title('corr(Wnew)')
subplot(224),imagesc(Wnew),caxis(max(abs(caxis))*[-1,1]); title('Wnew')

% figure,imexpb
% subplot(121),imagesc(W~=0),title(num2str(nnz(W)))
% subplot(122),imagesc(Wnew~=0),title(num2str(nnz(Wnew)))
figure,imexpb
subplot(121),imsupp(W)
subplot(122),imsupp(Wnew)
% return
%% sample stuffs
n = 500;
X = randn(n,p);
sig=1;
Y = X*Wnew + sig*randn(n,q);
Y2=X*W+sig*randn(n,q);
figure,imexpb
subplot(121),imcov(corr(Y2)), title('Unmixed weight vector')
subplot(122),imcov(corr(Y)),  title('Correlated output')

%% what if features have structured distribution?
X = tak_sample_AR1d(p,0.99,n);
Y = X*Wnew + sig*randn(n,q);
Y2=X*W+sig*randn(n,q);
figure,imexpb
subplot(121),imcov(corr(Y2)), title('Unmixed weight vector (AR features)')
subplot(122),imcov(corr(Y)),  title('Correlated output (AR features)')
