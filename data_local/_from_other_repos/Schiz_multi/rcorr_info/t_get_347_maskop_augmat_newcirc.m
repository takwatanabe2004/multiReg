%% t_get_347_maskop_augmat_newcirc.m 
%-------------------------------------------------------------------------%
% (02/12/2014)
% - same as t_get_347_maskop_augmat_ver2.m, but the scripts now assumes
%   that the "wrap-around-artifact" of the finite difference matrix takes 
%   place at the first row (previously, it was assumed to be in the last row).
% - found this approach to be far more convenient
%-------------------------------------------------------------------------%
% - half of the mask in the masking matrix created in t_get_347_maskop_augmat.m 
%   is redundance since the augmentation matrix fills up half of the vectors by
%   zeros...so this script replace these redundant masking entries by zeroes.
% - this does not change the penalty value of the graphnet or fused-lasso, but
%   during the ADMM operation, I've found some speedup can be attained by
%   having more 'zeroes' in the masking operator.
% create and save augmentation matrix and masking vector b.
% most of the script from t_j24_create_maskop_347yeo_CRIT
%%
clear all
purge

fsave=true;

outPath = [get_rootdir,'/data_local/rcorr_info/augmat_mask347newcirc.mat'];
outVars={'b', 'A', 'timeStamp',  'mFileName'};
%% this cell block is mostly identical to t_get_347_maskop_augmat.m
load graph_info347.mat adjmat coord
% randn('state',0)
d=coord.nsamp;
p=d*(d-1)/2;
NSIZE=coord.NSIZE;

N=coord.N;

%==================================================================================
% - Create augmentation matrix A (1st and 2nd level combined)
% - Also create array-converted connectome signal W, and see if it agrees with
%   reshape(A*w,[NSIZE NSIZE])
%==================================================================================
w=randn(p,1);

W=zeros(N,N);
A=sparse(N^2,p);

cnt=1;
for jj=1:d
    if mod(jj,50)==0; jj, end;
    for ii=jj+1:d
        ix=coord.rlex(ii);
        iy=coord.rlex(jj);
        W(ix,iy)=w(cnt);

        idx=sub2ind([N N],ix,iy);
        A(idx,cnt)=1;
        cnt=cnt+1;
    end
end

%=========================================================================%
% verify that the augmentation matrix does it's job, and also check
% if A'*A is identity
%=========================================================================%
if isequal( reshape(W,[NSIZE NSIZE]), reshape(A*w,[NSIZE NSIZE]) ) && isequal(speye(p),A'*A)
    disp('Good!!!')
else
    error('welp...screwed up again...now go enjoy debugging...')
end

% Apply brute force method create the neighborhood-graph (just as i did in t_create_6d_NN_adjmat.m)
%=========================================================================%
% difference matrix created using the adjacency matrix created from t_get_347_graphinfo.m
%=========================================================================%
C_brute=tak_adjmat2incmat(adjmat);
L_brute=C_brute'*C_brute;

%==================================================================================
% non-circulant difference matrix defined on the final augmented space
%==================================================================================
C_noncirc=tak_diffmat_newcirc([NSIZE NSIZE],0);
L_noncirc=C_noncirc'*C_noncirc;
% imedgel(C_noncirc,0),axis on

%==================================================================================
% circulant difference matrix defined on the final augmented space
%==================================================================================
C_circ=tak_diffmat_newcirc([NSIZE NSIZE],1);
L_circ=C_circ'*C_circ;

%==================================================================================
% some sanity check plots
%==================================================================================
% tspyl(L_brute)
% tspyl(L_noncirc)
% tspyl(L_circ)

idx_range=1:2000;
axis1=[idx_range(1),idx_range(end),idx_range(1),idx_range(end)];
% figure,hist(diag(L_circ)) % <- all coordinates have 12-neighbors here
% imcovvl(L_brute(idx_range,idx_range)),axis on,axis(axis1)
% imcovvl(L_noncirc(idx_range,idx_range)),axis on,axis(axis1)
% imcovvl(L_circ(idx_range,idx_range)),axis on,axis(axis1)
% imedger(L_brute(idx_range,idx_range),0),axis on,axis(axis1)
% imedger(L_noncirc(idx_range,idx_range),0),axis on,axis(axis1)
% imedger(L_circ(idx_range,idx_range),0),axis on,axis(axis1)

% obtain the desired difference and the graphnet & fused-lasso penalty
%==================================================================================
% the gold-standard penalty values from the brute-force difference matrix
%==================================================================================
% w=randn(p,1);
% w=ones(p,1);
load rcorr_design_censor X
w=X(1,:)';
W=reshape(A*w,[NSIZE,NSIZE]);
diff_brute=C_brute*w;
gnet_brute=norm(diff_brute,2)^2
flasso_brute=norm(diff_brute,1)

%==================================================================================
% create binary masking matrix to apply the circulant difference matrix
%==================================================================================
support_mask=(W~=0);

Bx1=circshift(support_mask,[+1  0  0  0  0  0])-support_mask;
By1=circshift(support_mask,[ 0 +1  0  0  0  0])-support_mask;
Bz1=circshift(support_mask,[ 0  0 +1  0  0  0])-support_mask;
Bx2=circshift(support_mask,[ 0  0  0 +1  0  0])-support_mask;
By2=circshift(support_mask,[ 0  0  0  0 +1  0])-support_mask;
Bz2=circshift(support_mask,[ 0  0  0  0  0 +1])-support_mask;

Bx1=tak_spdiag(Bx1(:)==0);
By1=tak_spdiag(By1(:)==0);
Bz1=tak_spdiag(Bz1(:)==0);
Bx2=tak_spdiag(Bx2(:)==0);
By2=tak_spdiag(By2(:)==0);
Bz2=tak_spdiag(Bz2(:)==0);

% blkdiag can be slow for large sparse matrices
NN=N^2;
Bsupp=...
[          Bx1, sparse(NN,NN), sparse(NN,NN), sparse(NN,NN), sparse(NN,NN), sparse(NN,NN); ...
 sparse(NN,NN),           By1, sparse(NN,NN), sparse(NN,NN), sparse(NN,NN), sparse(NN,NN); ...
 sparse(NN,NN), sparse(NN,NN),           Bz1, sparse(NN,NN), sparse(NN,NN), sparse(NN,NN); ...
 sparse(NN,NN), sparse(NN,NN), sparse(NN,NN),           Bx2, sparse(NN,NN), sparse(NN,NN); ...
 sparse(NN,NN), sparse(NN,NN), sparse(NN,NN), sparse(NN,NN),           By2, sparse(NN,NN); ...
 sparse(NN,NN), sparse(NN,NN), sparse(NN,NN), sparse(NN,NN), sparse(NN,NN),           Bz2];

%==================================================================================
% NOTE: unlike the 1d simulation in t_j23_1d_connectome_full_augmentation_CRIT.m,
%       the above support matrix must be composed with the binary circulant masker
%==================================================================================
Bcirc=tak_circmask_newcirc([NSIZE NSIZE]);
B=Bsupp*Bcirc;
b=logical(full(diag(B)));

if ~isdiagonal(B), error('meh...'), end;
tic
gnet_circ=norm(B*C_circ*W(:),2)^2;
flasso_circ=norm(B*C_circ*W(:),1);
toc

tic
gnet_circ2=norm(b.*(C_circ*W(:)),2)^2;
flasso_circ2=norm(b.*(C_circ*W(:)),1);
toc
isequal(gnet_circ,gnet_circ2)
isequal(flasso_circ,flasso_circ2)
Whos

[gnet_brute,gnet_circ]
[flasso_brute,flasso_circ]

err_gnet=abs(gnet_brute-gnet_circ)
err_flasso=abs(flasso_brute-flasso_circ)

if norm(gnet_brute-gnet_circ,1)<1e-8 && norm(flasso_brute-flasso_circ,1)<1e-8
% if isequal(gnet_brute,gnet_circ) && isequal(flasso_brute,flasso_circ)
    disp('I get same penalty value!!! SUCCESS! =)')
else
    error('what what???  how close were the two terms???')
end
% return
%% from here, the script diverges from the old one!
% remove 1's on the lexicographically upper-triangular part
mask=ones(N^2,1);

ix=zeros( N*(N-1)/2,1);
iy=zeros( N*(N-1)/2,1);
cnt=0;
for lex1=1:N; % 3d lex-ind
    for lex2=lex1+1:N % 3d lex-ind
        cnt=cnt+1;
        ix(cnt)=lex1;
        iy(cnt)=lex2;
    end
end
% convert from 2d-subs to 1d ind
idx=sub2ind([N N],ix,iy);
mask(idx)=0;
mask=tak_spdiag(mask);
% mask=reshape(mask,[NSIZE NSIZE]);
mask=blkdiag(mask,mask,mask,mask,mask,mask);
% mask=diag(mask);

B2=mask*Bsupp*Bcirc;
gnet_circ3=norm(B2*C_circ*W(:),2)^2;
flasso_circ3=norm(B2*C_circ*W(:),1);

%%%% sanity check %%%%
isequal(gnet_circ2,gnet_circ3)
isequal(flasso_circ2,flasso_circ3)


b=full(diag(B));
b2=full(diag(B2));

% the ratio of the 1's in the masking vector
sum(b)/length(b)
sum(b2)/length(b2)
%% save
timeStamp=tak_timestamp
mFileName=mfilename
clear b 
b=logical(b2);

if fsave
    save(outPath,outVars{:})
end
%%