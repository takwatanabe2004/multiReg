%% may09_get_maskop_cropped.m
% (05/08/2014)
%=========================================================================%
% - Find a way to obtain "masking matrix" B for the 3d case
% - code mostly from t_get_347_maskop_augmat_newcirc.m from my old project
%-------------------------------------------------------------------------%
% - this script w.r.t. "cropped FOV" 
%=========================================================================%
%%
clear
purge
disp('-------------------------------------------------------------------')


load([get_rootdir,'/data_local/mask_cropped.mat'],'nii_mask_crop','crop_range')
mask = nii_mask_crop.img;
[nx,ny,nz]=size(mask);

load([get_rootdir,'/data_local/diffmat_C.mat'],'C')
load([get_rootdir,'/data_local/A_matrix_cropped.mat'],'A')
% return
%%
L = C'*C;

NSIZE = [nx, ny, nz];
N = prod(NSIZE);
p=sum(mask(:));

%=========================================================================%
% circulant difference matrix defined on the final augmented space
%=========================================================================%
C_circ=tak_diffmat_newcirc(NSIZE,1);
L_circ=C_circ'*C_circ;
% tspy(L_circ)

%% obtain the desired difference and the graphnet & fused-lasso penalty
%=========================================================================%
% the gold-standard penalty values from the brute-force difference matrix
%=========================================================================%
randn('state',0); w=randn(p,1);
% w=ones(p,1);
% load sMRI_design_censor X
% w=X(1,:)';
W=reshape(A*w,NSIZE);
diff_brute=C*w;
gnet_brute=norm(diff_brute,2)^2;
flasso_brute=norm(diff_brute,1);

%=========================================================================%
% create binary masking matrix to apply the circulant difference matrix
%=========================================================================%
% support_mask=(W~=0);
support_mask=mask;
%%
% mask_nii2=mask_nii;
% mask_nii2.img=(W~=0);
% view_nii(mask_nii);
% view_nii(mask_nii2);
% return
%%
Bx=circshift(support_mask,[+1  0  0])-support_mask;
By=circshift(support_mask,[ 0 +1  0])-support_mask;
Bz=circshift(support_mask,[ 0  0 +1])-support_mask;
Bx=tak_spdiag(Bx(:)==0);
By=tak_spdiag(By(:)==0);
Bz=tak_spdiag(Bz(:)==0);

% blkdiag can be slow for large sparse matrices
Bsupp = [         Bx, sparse(N,N), sparse(N,N); ...
         sparse(N,N),          By, sparse(N,N); ...
         sparse(N,N), sparse(N,N),         Bz];
     
%==================================================================================
% NOTE: unlike the 1d simulation in t_j23_1d_connectome_full_augmentation_CRIT.m,
%       the above support matrix must be composed with the binary circulant masker
%==================================================================================
Bcirc=tak_circmask_newcirc(NSIZE);
B=Bsupp*Bcirc;
b=logical(full(diag(B)));

if ~isdiagonal(B), error('meh...'), end;
gnet_circ=norm(b.*(C_circ*W(:)),2)^2;
flasso_circ=norm(b.*(C_circ*W(:)),1);
% Whos

err_gnet=abs(gnet_brute-gnet_circ)
err_flasso=abs(flasso_brute-flasso_circ)
[size(A,1)/p, size(A,1)/1e6, length(b)/p, length(b)/1e6]