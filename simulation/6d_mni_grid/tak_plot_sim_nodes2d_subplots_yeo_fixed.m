function tak_plot_sim_nodes2d_subplots_yeo_fixed(H,nx,ny,idx_anom,idxNodes,...
    roiLabel,yeoColors,flagCoord)
% (06/30/2014)
%-------------------------------------------------------------------------%
% - fixed the error in tak_plot_sim_nodes2d_subplots_yeo.m and consequently
%   the error in june30_view_MNI_node_coverage_sliceBYslice_YeoColored.m...
%   the counter approach was not the right way....had to display the
%   lexico-index directly
% - also added flagCoord = {0,1,2,3}
%       * 0 = don't display xlabel & ylabel
%       * 1 = xlabel('x'), ylabel('y')
%       * 2 = xlabel('x'), ylabel('z')
%       * 3 = xlabel('y'), ylabel('z')
%-------------------------------------------------------------------------%
% - Display "counter" representing the lexicographic index of the nodes, 
%   as well as the color-code representing the yeoMembership of the nodes
% (run june30_view_MNI_node_coverage_sliceBYslice_YeoColored.m)
%-------------------------------------------------------------------------%
% (builds on top of tak_plot_sim_nodes2d_subplots3.m)
%=========================================================================%
% unlike version1, don't plot the "ghost nodes"
%-------------------------------------------------------------------------%
% - H = subplot handle
% - ultra makeshift version of tak_plot_sim_nodes2d.m for use in subplots
%   (where we need a specialized markersize......no intent to maintain this script)
%% some plot options
% get size of figure (wanna make the node-size proportional to it)
[tmp]=get(gcf,'position');
[tmp_subplot]=get(H,'position');
imsize_x=tmp_subplot(3)*tmp(3);
imsize_y=tmp_subplot(4)*tmp(4);
% keyboard
% return
imsize1 = min(imsize_x/nx,imsize_y/ny);
imsize2 = min(imsize_x,imsize_y);
% 
% msize=imsize/3000;
% mwidth=imsize/20000;
% 
% fsize=5
% lwidth = 5
%%
% msize=18;
% mwidth=4;

msize=imsize1/1.5;
mwidth=msize/8;

fsize=imsize2/20;
lwidth = imsize2/300;

% set background color
colorBack=[0.85 1 0.85];

axesOption={'XTick',[1,round(nx/2),nx],'YTick',[1,round(ny/2),ny],...
    'linewidth',lwidth,'Fontweight','b',...
    'fontsize',fsize,'TickLength',[0 0],'xaxislocation','bottom'};

% color_blue = [0.8 0.8 1];
% color_red = [1, 0.8 0.8];
color_blue = [0.6 0.6 1];
color_red = [1, 0.6 0.6];
markerOption={'o', 'MarkerEdgeColor','k','MarkerSize',msize,'linewidth',mwidth};
%% plot
imagesc(zeros(nx,ny)')

% set background color
colormap(colorBack) 

% display seeds in circular bubbles
hold on
axis('on','image','xy');
% [xx,yy]=ndgrid(1:nx,1:ny);
% plot(xx(:),yy(:),markerOption{:},'MarkerFaceColor', color_blue)

textOption={'fontsize',msize/2,'fontweight','b','HorizontalAlignment','center',...
            'VerticalAlignment','middle'};
%% plot nodes in brain coverage
%% if 3rd argument (anomalous node indices) is given, plot them in red
%---------------------------------------------------------------------% 
% (06/29/2014)
% here (x,y) indices of the anom-nodes are assumed to be given as input
%---------------------------------------------------------------------% 
if ~exist('flagCoord','var')
    flagCoord=0;
end

for i=1:size(idx_anom,1)
    ix=idx_anom(i,1);
    iy=idx_anom(i,2);
%         keyboard
    idx_yeo=roiLabel(idxNodes(i))+1; % +1 since roiLabel spans from 0~12
    plot(ix,iy,markerOption{:},'MarkerFaceColor', yeoColors(idx_yeo,:)    )
%         keyboard
%         str=['(',num2str(ix),',',num2str(iy),')'];
%         text(ix,iy,str,textOption{:},'fontsize',msize/2.5)
    text(ix,iy,num2str(idxNodes(i)),textOption{:})
    switch flagCoord
        case 1
            xlabel('y','fontsize',fsize*1.2,'fontweight','b')
            ylabel('z','fontsize',fsize*1.2,'fontweight','b')
        case 2
            xlabel('x','fontsize',fsize*1.2,'fontweight','b')
            ylabel('z','fontsize',fsize*1.2,'fontweight','b')
        case 3
            xlabel('x','fontsize',fsize*1.2,'fontweight','b')
            ylabel('y','fontsize',fsize*1.2,'fontweight','b')
    end
end
set(gca,axesOption{:})
axis xy
drawnow
% keyboard