%% july07_correlatedWeight_covapproach_patches
% (07/07/2014)
%=========================================================================%
% - same as july07_correlatedWeight_covapproach.m, but I made the weights
%   from each output have a patch structure.
%-------------------------------------------------------------------------%
% - The "covariance" approach to creating correlated output....
% - When the number of output grows, you won't get a sparse weight matrix...
%=========================================================================%
%%
clear all;
purge

% randn('state',0)
% rand('state',0)

%%
p = 10e3;
q = 10;

%% create (q x q) positive definite matrix
S = zeros(q,q);

%%%% for q = 10 %%%%%%%
% S(1:3,1:3)=randn(3);
% S(4:6,4:6)=randn(3);
% S(7:10,7:10)=randn(4);
% S=S*S';

S(1:3,1:3)=0.6;
S(4:6,4:6)=0.9;
S(7:10,7:10)=0.8;
S = S - diag(diag(S)) + eye(q); % <- replace diagonal entries with 1
%%%%%%%%%%%%%%%%%%%%%%%

% blk1 = 1:30;    nblk1=length(blk1);
% blk2 = 31:50;   nblk2=length(blk2);
% blk3 = 51:80;   nblk3=length(blk3);
% blk4 = 81:100;  bblk4=length(blk4);
% 
% % S(blk1,blk1)=randn(nblk1);
% S(blk1,blk1)=inv(tak_prec_ar1d(nblk1,0.95));
% % S(blk1,blk1)=.5*eye(nblk1) + 0.5;
% 
% % S(blk2,blk2)=randn(length(blk2));
% S(blk2,blk2)=inv(tak_prec_ar1d(length(blk2),0.9));
% 
% % S(blk3,blk3)=randn(length(blk3));
% S(blk3,blk3)=0.6*eye(length(blk3)) + 0.4;
% S(blk3,blk3)=inv(tak_prec_ar1d(length(blk3),0.9));
% 
% % S(blk4,blk4)=randn(length(blk4));
% % S(blk4,blk4)=0.8*eye(length(blk4)) + 0.2;
% S(blk4,blk4)=inv(tak_prec_ar1d(length(blk4),0.9));

% imcovvl(S)
L=chol(S,'lower');
% return
%% 
W = zeros(p,q);
for iq=1:q
    for ipatch = 1:5
        %=====================================================================%
        % create arbitrary patch for each outputs
        %=====================================================================%   
        patchLen = 80 * (5+randsample(5,1)); % patchLen \in K*{5,6,...,10}
        patchStart = randsample(p-patchLen,1); %
%         W(patchStart+1:patchStart+patchLen,iq) = tak_sample_signed_unif([5,10],1);
        tmp = zeros(p,1);
        tmp(patchStart+1:patchStart+patchLen) = tak_sample_signed_unif([5,10],1);
        W(:,iq) = W(:,iq) + tmp;
    end
end
Wnew = W*L;
figure,imexp
for ii=1:10
    subplot(3,4,ii),tplot(W(:,ii))
end
figure,imexp
for ii=1:10
    subplot(3,4,ii),tplot(Wnew(:,ii))
end

figure,imexp
subplot(221),imcov(corr(W),1), title('corr(W)')
subplot(223),imagesc(W),caxis(max(abs(caxis))*[-1,1]); title('W')
subplot(222),imcov(corr(Wnew),1), title('corr(Wnew)')
subplot(224),imagesc(Wnew),caxis(max(abs(caxis))*[-1,1]); title('Wnew')

figure,imexpb
subplot(121),imagesc(W~=0),title(num2str(nnz(W)))
subplot(122),imagesc(Wnew~=0),title(num2str(nnz(Wnew)))
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
