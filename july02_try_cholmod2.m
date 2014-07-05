clear
disp('----------------')

%-------------------------------------------------------------------------%
% p=50e3;
% C=tak_diffmat_1d(p,1);
% LAP  = C'*C + speye(p);
%-------------------------------------------------------------------------%
load('C:\Users\takanori\Documents\MATLAB\multiRegr\data_local\graphinfo\graph_info_Grid326.mat', 'C')
[~,p]=size(C);
LAP  = C'*C + 2*speye(p);


% [U,D,V]=svds(C,'econ');
% A = sparse(randn(5)+100*eye(5));
% cholmod2(A+A',randn(5,1))

tic
chol_L=lchol(LAP);
toc

% error('this method blows up')
% tic
% chol_LD=ldlchol(LAP);
% toc
figure,drawnow
%%
tic
niter=10;
for i=1:niter
chol_L\randn(p,1);
end
toc
toc/niter
%%
%-------------------------------------------------------------------------%
% method below enters swap! don't run
%-------------------------------------------------------------------------%
% tic
% error('DONT RUN THIS!! ENTERS SWAP!!!'); LL=chol(LAP')';
% toc