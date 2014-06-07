%% may30_spectr_rad_study.m
% (05/30/2014)
%=========================================================================%
% Studying simple spectral radiuses 
% (see my 2014 G-ADMM notes.pdf)...my handwritte notes
%=========================================================================%
%%
clear
purge
disp('***************************************************')
%% study identity matrices
p=40;
n=30;
% I = eye(p);
% disp('----- Concat horizontally ----')
% svds(I,1)
% svds([I,I],1)
% svds([I,I,I],1)
% disp('----- Concat vertically ----')
% svds(I,1)
% svds([I;I],1)
% svds([I;I;I],1)
%% study diff operators
% create row vector with single +1 and -1 entry at arbitrary point
alph=zeros(1,p);

% the random two entry locations (selection without replacement default)
ind=randsample(p,2);

alph(ind(1))=-1;
alph(ind(2))=+1;
% alph

svds(alph,1)
%%
%-------------------------------------------------------------------------%
% create concatenated alpha (vertical and horitzonaly)
%-------------------------------------------------------------------------%
alphHor=[];
alphVer=[];
for i=1:n
    alph=zeros(1,p);
    ind=randsample(p,2);
    alph(ind(1))=-1;
    alph(ind(2))=+1;
    alphHor=[alphHor,alph];
    alphVer=[alphVer;alph];
end
disp('----- Concat horizontally ----')
svds(alphHor,1)^2
disp('----- Concat vertically ----')
svds(alphVer,1)^2
%%
disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
[svds(alph,1)^2, n]
disp('---- repeat same diff horizontally ----')
alphHor=repmat(alph,[1,n]);
svds(alphHor,1)^2
disp('---- repeat same diff horizontally ----')
alphVer=repmat(alph,[n,1]);
svds(alphVer,1)^2
%%
disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
rank1=alphVer(1,:)'*alphVer(1,:);
rank2=alphVer(1:2,:)'*alphVer(1:2,:);
rank3=alphVer(1:3,:)'*alphVer(1:3,:);

[U1,V1]=eig(rank1);
V1=diag(V1)'
[U2,V2]=eig(rank2);
V2=diag(V2)'
[U3,V3]=eig(rank3);
V3=diag(V3)'

[U,V]=eig(alphVer'*alphVer);
V=diag(V)'
isequal(U1,U2,U3,U)