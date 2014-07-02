%% july02_loopless_corrDesignMatrix.m
% (07/02/2014)
%=========================================================================%
% - Given Z = (d x T x nsamp) data with d nodes, T time points, create
%   a loopless way to create design matrix X=(n x p), where p=nchoosek(d,2)
%   and X is a design matrix containing the correlations of Z
%=========================================================================%
%%
clear all;
purge

error('SIMPLE LOOP IS THE BEST...damn it...see SO link given in the script....')
%=========================================================================%
% sigh...so close...if only i could do A*A' for a 3-d array A in the 
% last cell block...
%-------------------------------------------------------------------------%
% ref: 
% - http://stackoverflow.com/questions/6580656/matlab-how-to-vector-multiply-two-arrays-of-matrices
%=========================================================================%

nsamp=30;
d=20;
T=50;

p=nchoosek(d,2);

Z = randn(d,T,nsamp);

%% loop way
X1 = zeros(nsamp,p);
for i=1:nsamp
    tmp = corr(Z(:,:,i)');
    X1(i,:) = tak_dvec(tmp)';
end
%% loopless way?

%=========================================================================%
% try for a single slice
%=========================================================================%
purge
i=1;
zslice = Z(:,:,i);

tmp=bsxfun(@minus, zslice, mean(zslice,2)); % zero-mean
tmp2=bsxfun(@rdivide,tmp, sqrt(T*var(zslice,1,2))     ); % normalize
% tmp2=bsxfun(@times,tmp, 1./  sqrt( sum(zslice.^2,2)  )     ); % normalize
mean(tmp,2)
var(tmp2,1,2)

norm(corr(zslice') - tmp2*tmp2')
% imcovvl(corr(zslice'))
% imcovvl(tmp2*tmp2')
%% loopless way: alll samples at once?

%-------------------------------------------------------------------------%
% zero mean
%-------------------------------------------------------------------------%
tmp1=mean(Z,2);
TMP1=bsxfun(@minus, Z, tmp1); % zero-mean
% tplottl(tak_vec(mean(TMP1,2)))

%-------------------------------------------------------------------------%
% normalize
%-------------------------------------------------------------------------%
tmp2 = var(TMP1,1,2);
TMP2=bsxfun(@rdivide,TMP1, sqrt(T*tmp2)     ); % normalize
% tplottl(tak_vec(var(TMP2,1,2)))
% imcovvl(corr(zslice'))
% imconnl(X1(1,:))
% imcovvl(corr(TMP2(:,:,1)'))
% imcovvl(TMP2(:,:,1)*TMP2(:,:,1)')

A = TMP2;
%-------------------------------------------------------------------------%
At = permute(A,[2 1 3]);
C = zeros(d,d,nsamp);
% C(1,1,:) = A(1,1,:).*At(1,1,:) + A(1,2,:).*At(2,1,:);
% C(1,2,:) = A(1,1,:).*At(1,2,:) + A(1,2,:).*At(2,2,:);
% C(2,1,:) = A(2,1,:).*At(1,1,:) + A(2,2,:).*At(2,1,:);
% C(2,2,:) = A(2,1,:).*At(1,2,:) + A(2,2,:).*At(2,2,:);
C(1,1,:) = A(1,1,:).*At(1,1,:) + A(1,2,:).*At(2,1,:);
C(1,2,:) = A(1,1,:).*At(1,2,:) + A(1,2,:).*At(2,2,:);
C(2,1,:) = A(2,1,:).*At(1,1,:) + A(2,2,:).*At(2,1,:);
C(2,2,:) = A(2,1,:).*At(1,2,:) + A(2,2,:).*At(2,2,:);
imcovvl(A(:,:,i)*At(:,:,i))
imcovvl(C(:,:,i))
% arrayfun(At,

% A2=TMP2.^2;
% A3=sum(A2,2)
% TMP3 = bsxfun(@times, TMP2, TMP2')
