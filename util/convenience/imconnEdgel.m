function imconnEdgel(w)
% imconnEdge(Im,flag)
%=========================================================================%
% - reshape vectorized connectome into symmetric matrix
%=========================================================================%
% (07/01/2014)
%%
figure,imexpl
wmat = tak_dvecinv(w,0);
imedge(wmat)

titleOption = {'Interpreter','none','fontweight','b','fontsize',16};
[nnz_upper, sp_level] = tak_nnz_lower(wmat);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_upper,sp_level*100);
title([inputname(1), ' --- ', titleStr],titleOption{:})
drawnow
