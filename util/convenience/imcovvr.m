function imcovvr(Im,flag)
% (07/05/2014) - added "flag" to optionally "symmetrize" the colormap scale
%                (default: off)
%%
if nargin == 1
    flag = 0;
end
figure,imexpl
imcov(Im,flag)

%| http://www.mathworks.com/support/solutions/en/data/1-16FLM/?solution=1-16FLM
title(inputname(1),'Interpreter','none')
drawnow