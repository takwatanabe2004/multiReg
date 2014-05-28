function imtriag(im)
% imtriag(im)
%|------------------------------------------------------------------------------|%
%|  Generate trinary representation of the result
%| (red = positive, blue = negative, white = zero)
%|------------------------------------------------------------------------------|%
%| (12/22/2013)
%%
CMAP=[0 0 1; 1 1 1; 1 0 0];

edgeGraph=im;
edgeGraph(edgeGraph>0)=+1; % positive = red
edgeGraph(edgeGraph<0)=-1; % negative = blue
edgeGraph(edgeGraph==0)=0;

imagesc(edgeGraph),colormap(CMAP)
axis('on','image')
[nnz_upper, sp_level] = tak_nnz_upper(im);
titleStr = sprintf('%g nonzeroes (%g%% sparsity level)', nnz_upper,sp_level*100);
title([inputname(1), ' --- ', titleStr],'Interpreter','none')
drawnow
% impixelinfo

% cbarOption={'fontsize',22','fontweight','b','ytick',[-.66,0,.66],...
%     'YTickLabel',{' <0',' =0',' >0'},'TickLength',[0 0]};
% colorbar(cbarOption{:})