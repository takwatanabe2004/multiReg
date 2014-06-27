clear
purge
p=100;
n=100;
q=10;
X=randn(n,p);
XX = kron(speye(q),X);
% figure,imexpl,imagesc(X)
% figure,imexpl,imagesc(XX),drawnow

svds(X,1)
svds(XX,1)




