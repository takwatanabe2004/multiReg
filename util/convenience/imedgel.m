function imedgel(im, zeroDiag)
% imedgel(im,zeroDiag)
%|------------------------------------------------------------------------------|%
%|  Convenient version of imedge.m...
%|------------------------------------------------------------------------------|%
%| 09/23/2012
%| 06/18/2013 -> impixelinfo added
%%
if(~exist('zeroDiag','var')||isempty(zeroDiag)), zeroDiag = true; end
figure,imexpl
%| adjacency matrix
adjMat = im~=0;

if zeroDiag == true
    adjMat(logical(eye(size(adjMat)))) = false;
end

imagesc(adjMat~=0)
axis('off','image')
colormap(flipud(gray))

[nnz_upper, sp_level] = tak_nnz_upper(im);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_upper,sp_level*100);
title([inputname(1), ' --- ', titleStr],'Interpreter','none')
drawnow
impixelinfo