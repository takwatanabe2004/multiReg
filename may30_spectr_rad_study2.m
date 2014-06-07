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
p=5;
nrep=3;
I = eye(p);
disp('----- Concat horizontally ----')
Ihor2=[I,I];
Ihor3=[I,I,I];
svds(Ihor2,1)
svds(Ihor3,1)
disp('----- Concat vertically ----')
Iver2=[I,I]';
Iver3=[I,I,I]';
svds(Iver2,1)
svds(Iver3,1)
return
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
for i=1:nrep
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
% disp('---- repeat same diff horizontally ----')
% alphHor=repmat(alph,[1,nrep]);
% svds(alphHor,1)
% disp('---- repeat same diff horizontally ----')
% alphVer=repmat(alph,[nrep,1]);
% svds(alphVer,1)