function imcov(Im,flag)
% (07/05/2014) - added "flag" to optionally "symmetrize" the colormap scale
%                (default: off)
%%
if nargin == 1
    flag = 0;
end

imagesc(Im)

if flag
    % symmetrize colormap
    tmp=caxis;
    tmp=max( abs(tmp) );
    caxis(tmp*[-1,1])
end

%| http://www.mathworks.com/support/solutions/en/data/1-16FLM/?solution=1-16FLM
title(inputname(1),'Interpreter','none')
axis('on','image');
% colorbar
colorbar('location','southoutside')
impixelinfo
drawnow