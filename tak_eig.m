function [U,V]=tak_eig(K)
% (05/30/2014)
% cuz matlab sorts eigvals in ascending order...
[U,V]=eig(K);
[V,idx]=sort(diag(V),'descend');
U=U(:,idx);
