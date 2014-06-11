% (06/10/2014)
%=========================================================================%
% study feature matrix
%=========================================================================%
clear
purge

grid='Grid326'; % {'Grid326','Grid1068','WashU'}


%%
% visualization of the result
purge
cbarOption={'fontsize',22','fontweight','b','ytick',[-.66,0,.66],...
    'YTickLabel',{' <0',' =0',' >0'},'TickLength',[0 0]};
textOption1={'fontweight','b','fontsize',9};
lineOption = {'color','k','linewidth',0.5};
lwidth_deg=2.5;
%%
load([get_rootdir,'/data_local/designMatrix_FC_',grid,'.mat'])
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'])

% circularly shift 1 indices (so "unlabeled" is at the final label index)
roiLabel=roiLabel-1;
roiLabel(roiLabel==-1)=12;
yeoLabels=circshift(yeoLabels,-1);
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);

xmat1 = tak_dvecinv(X(21,:),0);
xmat2 = tak_dvecinv(X(22,:),0);

xmat1_srt=xmat1(idxsrt,idxsrt);
xmat2_srt=xmat2(idxsrt,idxsrt);
% sum(isnan(X(:)))
% imcovvl()
% imcovvl(tak_dvecinv(X(2,:),0))
% 
%%
cmap=[-1,1];
imcovvl(xmat1_srt),axis off,colorbar('location','northoutside'),caxis(cmap)
tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)

imcovvl(xmat2_srt),axis off,colorbar('location','northoutside'),caxis(cmap)
tak_local_linegroups5(gcf,labelCount,textOption1,yeoLabels,lineOption)