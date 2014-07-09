function imsuppl(W)
% imsupp(W)
%=========================================================================%
% - Comments
%=========================================================================%
% (07/07/2014)
%%
figure,imexpl
imagesc(W~=0)
% axis('off','image') 
% axis('image') 
colormap(flipud(gray))


spLevel = nnz(W)/numel(W);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz(W),spLevel*100);
title([inputname(1), ' --- ', titleStr],'Interpreter','none')
drawnow
% impixelinfo