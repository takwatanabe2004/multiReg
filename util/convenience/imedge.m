function imedge(im, zeroDiag)
% imedge(im,zeroDiag)
%|------------------------------------------------------------------------------|%
%|  Display the support of the image, ie, the non-zero entries.
%|  Use black & white display, with black being the nonzero components.
%|
%|  If 'zeroDiag' = true, the diagonal of the adjacency matrix will be zeroed out.
%|     (default = true)
%|------------------------------------------------------------------------------|%
% (06/27/2014) -> modified tak_nnz_upper to tak_nnz_lower
%| 09/07/2012
%| 06/18/2013 -> impixelinfo added
%%
if(~exist('zeroDiag','var')||isempty(zeroDiag)), zeroDiag = true; end

%| adjacency matrix
adjMat = im~=0;

if zeroDiag == true
    adjMat(logical(eye(size(adjMat)))) = false;
end

imagesc(adjMat~=0)
% axis('off','image') 
axis('image') 
colormap(flipud(gray))

[nnz_upper, sp_level] = tak_nnz_lower(im);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_upper,sp_level*100);
title([inputname(1), ' --- ', titleStr],'Interpreter','none')
drawnow
% impixelinfo