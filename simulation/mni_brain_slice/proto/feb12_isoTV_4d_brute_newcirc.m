%% feb12_isoTV_4d_brute_newcirc.m
% (02/12/2014)
%=========================================================================%
% Compute the isotropic TV value of the 4-D data by brute force
% - use the same node orientation as the single task simulation
% - use the new circulant matrix condition that assumes the first row
%   contains the "wrap-around effect."
%=========================================================================%
%%
clear
purge
load('graph_info347_2d.mat','coord','adjmat')
load('augmat_mask347newcirc_2d.mat', 'A','b')
% randn('state',0)
[ptil,p]=size(A);
w=randn(p,1);
wtil=A*w;
C=tak_adjmat2incmat(adjmat);
Cw=C*w;
Aw=A*w;
%% get operators of interest (gradient, masking, circulant gradient, etc)
%-------------------------------------------------------------------------%
% collect the gradient operator for each connexel location
%-------------------------------------------------------------------------%
D=cell(p,1); 
for i=1:p
    idx=(C(:,i)==1);
    D{i}=C(idx,:);
end

%-------------------------------------------------------------------------%
% circulant gradient operator at the matrix level
%-------------------------------------------------------------------------%
[C_circ, Cx1_circ, Cy1_circ, Cx2_circ, Cy2_circ] = ...
    tak_diffmat_4d_newcirc([coord.NSIZE,coord.NSIZE],1); % <- C'*C has circulant structure!
% figure,imexp
% subplot(2,3,[1,4]),tspy(C_circ)
% subplot(232),tspy(Cx1_circ)
% subplot(233),tspy(Cx2_circ)
% subplot(235),tspy(Cy1_circ)
% subplot(236),tspy(Cy2_circ)
% return
%-------------------------------------------------------------------------%
% laplacian matrices
%-------------------------------------------------------------------------%
Lap = C'*C;
Lap_circ = C_circ'*C_circ;

%-------------------------------------------------------------------------%
% Create the gradient operator for each connexel of the augmented connectome
%-------------------------------------------------------------------------%
Dtil=cell(ptil,1); 
btil=cell(ptil,1);
for i=1:ptil
    idx=(C_circ(:,i)==1);
    Dtil{i}=C_circ(idx,:);
    btil{i}=b(idx);
%     [i,sum(idx)]
end
% return
%%
disp('---- Sanity check on the fused lasso values ----')
% matrix-vector form
flasso_matvec=norm(C*w,1);

% summation form
flasso_summ1=0;
for i=1:p
    if ~isempty(D{i})
        flasso_summ1 = flasso_summ1 + norm(D{i}*w,1);
    end
end
flasso_matvec
flasso_summ1

%-------------------------------------------------------------------------%
% repeat the summation form above without using D{i}'s
%-------------------------------------------------------------------------%
% summation form
flasso_summ2=0;
for i=1:p
    idx=find(C(:,i)==1);
    if ~isempty(idx)
        flasso_summ2 = flasso_summ2 + norm(Cw(idx),1);
    end
end
flasso_summ2
if isequal(flasso_summ1,flasso_summ2)
    disp('Good...the fused lasso value without creating Di''s agreed')
else
    error('crap...')
end
% return

%-------------------------------------------------------------------------%
% Now for the AUGMENTED connectome
%-------------------------------------------------------------------------%
% matrix-vector form on the augmented connectome
flasso_matvec_circ = norm(b.*(C_circ*Aw),1)

% "summation form" on the augmented connectome
flasso_summ_circ1=0;
for i=1:ptil
    flasso_summ_circ1 = flasso_summ_circ1 + ...
        norm(btil{i}.*(Dtil{i}*wtil),1);
end
flasso_summ_circ1

%-------------------------------------------------------------------------%
% repeat the summation form without using btil{i}'s and Dtil{i}'s
%-------------------------------------------------------------------------%
CAw=reshape(C_circ*Aw,[4,ptil]);
bmat=reshape(b, [4,ptil]);
flasso_summ_circ2=0;
for i=1:ptil
    flasso_summ_circ2 = flasso_summ_circ2 + ...
        norm(bmat(:,i).*CAw(:,i),1);
end
err=abs(flasso_summ_circ1 - flasso_summ_circ2)
if err < 1e-10
    disp('Good...the fused lasso value without creating Dtil{i}''s agreed')
else
    error('crap...')
end
% return

%=========================================================================%
% evaluate accuracy/error
%=========================================================================%
err12=abs(flasso_matvec-flasso_summ1);
err13=abs(flasso_matvec-flasso_matvec_circ);
err14=abs(flasso_matvec-flasso_summ_circ1);
err23=abs(flasso_summ1-flasso_matvec_circ);
err24=abs(flasso_summ1-flasso_summ_circ1);
err34=abs(flasso_matvec_circ-flasso_summ_circ1);
err=[err12,err13,err14,err23,err24,err34]
if sum(err) < 1e-8
    disp('Success with fused lasso values!!!')
else
    error('welp...go debug')
end
% return
%% repeat above for graphnet
disp('---- Sanity check on the GraphNet values ----')
% matrix-vector form
graphNet_matvec=norm(C*w,2)^2

% summation form
graphNet_summ1=0;
for i=1:p
    if ~isempty(D{i})
        graphNet_summ1 = graphNet_summ1 + norm(D{i}*w,2)^2;
    end
end
graphNet_summ1

%-------------------------------------------------------------------------%
% repeat the summation form above without using D{i}'s
%-------------------------------------------------------------------------%
graphNet_summ2=0;
for i=1:p
    idx=find(C(:,i)==1);
    if ~isempty(idx)
        graphNet_summ2 = graphNet_summ2 + norm(Cw(idx),2)^2;
    end
end
graphNet_summ2
if isequal(graphNet_summ1,graphNet_summ2)
    disp('Good...the GraphNet value without creating Di''s agreed')
else
    error('crap...')
end
% return

% matrix-vector form on the augmented connectome
graphNet_matvec_circ = norm(b.*(C_circ*A*w),2)^2

% "summation form" on the augmented connectome
graphNet_summ_circ1=0;
for i=1:ptil
    graphNet_summ_circ1 = graphNet_summ_circ1 + ...
        norm(btil{i}.*(Dtil{i}*wtil),2)^2;
end
graphNet_summ_circ1

%-------------------------------------------------------------------------%
% repeat the summation form without using btil{i}'s and Dtil{i}'s
%-------------------------------------------------------------------------%
CAw=reshape(C_circ*Aw,[4,ptil]);
bmat=reshape(b, [4,ptil]);
graphNet_summ_circ2=0;
for i=1:ptil
    graphNet_summ_circ2 = graphNet_summ_circ2 + ...
        norm(bmat(:,i).*CAw(:,i),2)^2;
end
err=abs(graphNet_summ_circ1 - graphNet_summ_circ2)
if err < 1e-10
    disp('Good...the fused lasso value without creating Dtil{i}''s agreed')
else
    error('crap...')
end
% return

err12=abs(graphNet_matvec-graphNet_summ1);
err13=abs(graphNet_matvec-graphNet_matvec_circ);
err14=abs(graphNet_matvec-graphNet_summ_circ1);
err23=abs(graphNet_summ1-graphNet_matvec_circ);
err24=abs(graphNet_summ1-graphNet_summ_circ1);
err34=abs(graphNet_matvec_circ-graphNet_summ_circ1);
err=[err12,err13,err14,err23,err24,err34]
if sum(err) < 1e-9
    disp('Success with GraphNet values!!!')
else
    error('welp...go debug')
end
%% compute isotropic TV values using summation & summ+mask method
disp('---- Sanity check on isotropic TV values ----')
%-------------------------------------------------------------------------%
% compute fused lasso via the "summation form"
%-------------------------------------------------------------------------%
isoTV1=0;
for i=1:p
    if ~isempty(D{i})
        isoTV1 = isoTV1 + norm(D{i}*w,2);
    end
end
isoTV1

%-------------------------------------------------------------------------%
% repeat the summation form above without using D{i}'s
%-------------------------------------------------------------------------%
Cw=C*w;
isoTV2=0;
for i=1:p
    idx=find(C(:,i)==1);
    if ~isempty(idx)
        isoTV2 = isoTV2 + norm(Cw(idx),2);
    end
end
isoTV2

if isequal(isoTV1,isoTV2)
    disp('Good...the isoTV value without creating Di''s agreed')
else
    error('crap...')
end
% return

%-------------------------------------------------------------------------%
% compute isotropic TV with the augmented signal
%-------------------------------------------------------------------------%
isoTV_circ1=0;
for i=1:ptil
    isoTV_circ1 = isoTV_circ1 + norm(btil{i}.*(Dtil{i}*A*w),2);
end
isoTV_circ1

%-------------------------------------------------------------------------%
% repeat the summation form without using btil{i}'s and Dtil{i}'s
%-------------------------------------------------------------------------%
isoTV_circ2=0;
CAw=C_circ*Aw;
dx1=CAw(1       :   ptil)';
dy1=CAw(1+  ptil: 2*ptil)';
dx2=CAw(1+2*ptil: 3*ptil)';
dy2=CAw(1+3*ptil: 4*ptil)';
GRADIENT = [dx1;dy1;dx2;dy2]; 
bmat=[b(1       :   ptil)'; b(1+  ptil: 2*ptil)'; ...
      b(1+2*ptil: 3*ptil)'; b(1+3*ptil: 4*ptil)'];
for i=1:ptil
    isoTV_circ2 = isoTV_circ2 + ...
        norm(bmat(:,i).*GRADIENT(:,i),2);
end
isoTV_circ2
err=abs(isoTV_circ1 - isoTV_circ2)

if sum(err) < 1e-10
    disp('Good...the isoTV value without creating Dtil{i}''s agreed')
else
    error('welp...go debug')
end

err1=abs(isoTV1-isoTV_circ1)
err2=abs(isoTV1 - isoTV2)
err3=abs(isoTV_circ1 - isoTV_circ2)
if sum(err1) < 1e-10
    disp('Success with isoTV values!!!')
else
    error('welp...go debug')
end
%%
%-------------------------------------------------------------------------%
% sanity check on gradient from diffmat and circshift
%-------------------------------------------------------------------------%
% figure,imexp
% subplot(121),imcov(reshape(dx,[coord.nx,coord.ny])')
% subplot(122),imcov(reshape(dy,[coord.nx,coord.ny])')

Wtil=reshape(Aw,[coord.NSIZE,coord.NSIZE]);
dx1_ver2 = Wtil - circshift(Wtil,[+1  0  0  0]);
dy1_ver2 = Wtil - circshift(Wtil,[ 0 +1  0  0]);
dx2_ver2 = Wtil - circshift(Wtil,[ 0  0 +1  0]);
dy2_ver2 = Wtil - circshift(Wtil,[ 0  0  0 +1]);
if isequal(dx1(:),dx1_ver2(:)) && isequal(dy1(:),dy1_ver2(:)) && ...
   isequal(dx2(:),dx2_ver2(:)) && isequal(dy2(:),dy2_ver2(:))
    disp('Finite diffmat and circshift gives same gradient maps')
else
    error('welp...debugging time...')
end
%% check if summing through the difference operator returns the BCCB lap mat
%=========================================================================%
% -> this takes a long time, but I confimed this is true.
%=========================================================================%
% disp('------- Laplacian summation check --------')
% Lap_circ_summ=zeros(ptil,ptil);
% tic
% for i=1:ptil
%     Lap_circ_summ = Lap_circ_summ + Dtil{i}'*Dtil{i};
% end
% if isequal(Lap_circ_summ, Lap_circ)
%     disp('success!')
% end