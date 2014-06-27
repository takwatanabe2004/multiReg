function imconn(w,arg2,arg3)
% imConn(Im,arg2,arg3)
%=========================================================================%
% - reshape vectorized connectome into symmetric matrix
%=========================================================================%
% (06/27/2014)
%%
wmat = tak_dvecinv(w,0);
if nargin==1
    imcov(wmat)
elseif nargin==2
    imcov(arg2)
else
    imcov(Im,[arg3(1) arg3(2)])
end
title(inputname(1),'Interpreter','none')
drawnow