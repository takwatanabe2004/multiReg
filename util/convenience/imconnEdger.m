function imconnEdger(w)
% imconnEdge(Im,flag)
%=========================================================================%
% - reshape vectorized connectome into symmetric matrix
%=========================================================================%
% (07/01/2014)
%%
figure,imexpr
wmat = tak_dvecinv(w,0);
imedge(wmat)

titleOption = {'Interpreter','none','fontweight','b','fontsize',16};
[nnz_lower, sp_level] = tak_nnz_lower(wmat);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_lower,sp_level*100);
title([inputname(1), ' --- ', titleStr],titleOption{:})
drawnow
