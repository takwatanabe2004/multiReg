%% may30_spectr_rad_study3.m
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
p=50;
nadj=4;
if nadj>=p
    error('nadj must be less than p')
end
% I = eye(p);
% disp('----- Concat horizontally ----')
% svds(I,1)
% svds([I,I],1)
% svds([I,I,I],1)
% disp('----- Concat vertically ----')
% svds(I,1)
% svds([I;I],1)
% svds([I;I;I],1)
%% study ADJACENT diff operators
A=zeros(nadj,p);
for i=1:nadj
    A(i,i)=-1;
    A(i,i+1)=+1;
end

% return
% create row vector with single +1 and -1 entry at arbitrary point
alph=A(1,:);
svds(alph,1)^2

svds(A,1)^2
AtA=A'*A;
[U,V]=tak_eig(AtA);
V'
figure,imexptr
imcov(U)
%%
C=full(tak_diffmat_1d(p,1));
[UU,VV]=tak_eig(C'*C);
VV'
figure,imexptr
imcov(UU)
%%
% C=full(tak_diffmat_1d(p));
% [UU,VV]=tak_eig(C'*C);
% VV'
% figure,imexptr
% imcov(UU)
%%
% disp('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
% rank1=alphVer(1,:)'*alphVer(1,:);
% rank2=alphVer(1:2,:)'*alphVer(1:2,:);
% rank3=alphVer(1:3,:)'*alphVer(1:3,:);
% 
% [U1,V1]=eig(rank1);
% V1=diag(V1)'
% [U2,V2]=eig(rank2);
% V2=diag(V2)'
% [U3,V3]=eig(rank3);
% V3=diag(V3)'
% 
% [U,V]=eig(alphVer'*alphVer);
% V=diag(V)'
% isequal(U1,U2,U3,U)