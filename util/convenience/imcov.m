% use this to have the axis imaging for covariance, default colormap
function imcov(Im,arg2,arg3)
if nargin == 1
    imagesc(Im)
elseif nargin == 2 
    if isa(arg2,'char') || (numel(arg2) > 2)
        imagesc(Im)
        colormap(arg2)
    else
        imagesc(Im,[arg2(1) arg2(2)])
        colormap jet
    end
else
    imagesc(Im,[arg3(1) arg3(2)])
    colormap(arg2)
end
%| http://www.mathworks.com/support/solutions/en/data/1-16FLM/?solution=1-16FLM
title(inputname(1),'Interpreter','none')
axis('on','image');
colorbar
impixelinfo
drawnow