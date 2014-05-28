% use this to have the axis imaging for covariance, default colormap
function imcovvl(Im,arg2,arg3)
figure
imexpl
if nargin == 1
    imcov(Im)
elseif nargin == 2
    imcov(Im, arg2)
else
    imcov(Im, arg2, arg3)
end
%| http://www.mathworks.com/support/solutions/en/data/1-16FLM/?solution=1-16FLM
title(inputname(1),'Interpreter','none')
drawnow