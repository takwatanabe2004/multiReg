function imconnEdge(w)
% imconnEdge(Im,flag)
%=========================================================================%
% - reshape vectorized connectome into symmetric matrix
%=========================================================================%
% (06/29/2014)
%%
wmat = tak_dvecinv(w,0);
imedge(wmat)

[nnz_upper, sp_level] = tak_nnz_lower(wmat);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_upper,sp_level*100);
title([inputname(1), ' --- ', titleStr],'Interpreter','none')
drawnow
