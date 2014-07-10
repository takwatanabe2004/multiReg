%% july09_corrWeight_covMethod.m
% (07/09/2014)
%=========================================================================%
% - Create sparse multi-output weight vector with correlated columns.
% - Use "covariance" approach.
%=========================================================================%
%%
clear all;
purge
error('DITCHED THIS METHOD!  THE CHOLESKY OR EIG BASED MATRIX SQUARED ROOT DOES A POOR JOB SMOOTHING FEATURES ACROSS CORRELATED OUTPUTS!!!')
randn('state',0)
rand('state',0)

p = 500;
q = 30;

% number of output block structure
nBlocks = 4;
%% create (q x q) PD matrix assiginig the crr. struct of the weight vectors
S = zeros(q,q);

idxCluster = cell(nBlocks,1);
% corrBlkList = [0.9,-.4,0.8,-0.3];
corrBlkList = 0.8*ones(nBlocks,1);
%%%% for numBlocks clusters of correlated features %%%%%%
offset = 0;
for iBlock = 1:nBlocks
    if iBlock == nBlocks
        idxCluster{iBlock} = offset+1:q;
    else
        idxCluster{iBlock} = offset+1:offset+round(q/nBlocks);
        offset = idxCluster{iBlock}(end);
    end
    
    idx = idxCluster{iBlock};
%     S(idx,idx)=corrBlkList(iBlock);
    S(idx,idx) = inv(tak_prec_ar1d( length(idx), sign(1)*0.99));
end
S = S - diag(diag(S)) + eye(q); % <- replace diagonal entries with 1

[~,flagPD]=chol(S);
%%%%%%%% Add diagonals until matrix is positive definite %%%%%%%%
cnt=1;
while(flagPD~=0)
    S = S + 0.05*eye(q);
    [~,flagPD]=chol(S);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=chol(S,'lower');
% [U,D]=eig(S);
% L = U*sqrt(D);


figure,imexpb
subplot(131),imcov(S)
subplot(132),imcov(L)
subplot(133),imcov(L*L')
% [~,spectrum]=tak_eig(S); tstem(spectrum) % <- check PD condition of S
% return
%%%%%%%%%%%%%%%%%%%%%%%
%% 
W = zeros(p,q);
nnz_col = round(p/10); % <- # nnz in each columns
for iq=1:q
    idx = randsample(p,nnz_col);
    W(idx,iq) = tak_sample_signed_unif([5,10],nnz_col);
end
Wmixed = W*L;
figure,imexp
for ii=1:min(q,10)
    subplot(3,4,ii),tplot(W(:,ii))
end
figure,imexp
for ii=1:min(q,10)
    subplot(3,4,ii),tplot(Wmixed(:,ii))
end
% return
figure,imexp
subplot(221),imcov(corr(W),1), title('corr(W)'), colorbar off, colorbar Eastoutside
subplot(223),imagesc(W),caxis(max(abs(caxis))*[-1,1]); title('W')
subplot(222),imcov(corr(Wmixed),1), title('corr(Wnew)'), colorbar off, colorbar Eastoutside
subplot(224),imagesc(Wmixed),caxis(max(abs(caxis))*[-1,1]); title('Wnew')

% figure,imexpb
% subplot(121),imagesc(W~=0),title(num2str(nnz(W)))
% subplot(122),imagesc(Wnew~=0),title(num2str(nnz(Wnew)))
figure,imexpb
subplot(121),imsupp(W)
subplot(122),imsupp(Wmixed)
%%
% purge
figure,imexp
for i=1:nBlocks
    subplot(4,1,i),tstem(Wmixed(:,idxCluster{i}(1:3)))
end
% return
%% sample stuffs
n = 500;
X = randn(n,p);
sig=1;
Y = X*Wmixed + sig*randn(n,q);
Y2=X*W+sig*randn(n,q);
figure,imexpb
subplot(121),imcov(corr(Y2)), title('Unmixed weight vector')
subplot(122),imcov(corr(Y)),  title('Correlated output')

%% what if features have structured distribution?
X = tak_sample_AR1d(p,0.99,n);
Y = X*Wmixed + sig*randn(n,q);
Y2=X*W+sig*randn(n,q);
figure,imexpb
subplot(121),imcov(corr(Y2)), title('Unmixed weight vector (AR features)')
subplot(122),imcov(corr(Y)),  title('Correlated output (AR features)')
