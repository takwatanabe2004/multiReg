function imconnr(w,flag)
% imconnr(w,flag)
%=========================================================================%
% - reshape vectorized connectome into symmetric matrix
% - (06/29/2014) flag will "symmetrize" the colormap scale (default: off)
%=========================================================================%
% (06/27/2014)
%%
figure,imexpr
if nargin == 1
    flag = 0;
end
wmat = tak_dvecinv(w,0);
imcov(wmat)
if flag
    tmp=caxis;
    tmp=max( abs(tmp) );
    caxis(tmp*[-1,1])
end
title(inputname(1),'Interpreter','none')
drawnow