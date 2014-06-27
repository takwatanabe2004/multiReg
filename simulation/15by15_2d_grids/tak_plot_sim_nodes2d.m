function tak_plot_sim_nodes2d(nx,ny,idx_anom)
% tak_function
%=========================================================================%
% - Given (nx by ny) node orientation as input, plot the nodes in 2d space.
% - 3rd argument (idx_anom) are anomalous nodes, plotted in red
%=========================================================================%
% (06/26/2014)
%% some plot options
% get size of figure (wanna make the node-size proportional to it)
[tmp]=get(gcf,'position');
imsize_x=tmp(3)
imsize_y=tmp(4)
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

msize=imsize1/2;
mwidth=msize/8;

fsize=imsize2/50;
lwidth = imsize2/300;

% set background color
colorBack=[0.85 1 0.85];

axesOption={'XTick',[1:5:nx-1,nx],'YTick',[1:5:ny-1,ny],...
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
axis('on','image');
% [xx,yy]=ndgrid(1:nx,1:ny);
% plot(xx(:),yy(:),markerOption{:},'MarkerFaceColor', color_blue)

textOption={'fontsize',msize/2,'fontweight','b','HorizontalAlignment','center',...
            'VerticalAlignment','middle'};
cnt=1;
for j=1:ny
    for i=1:nx
        plot(i,j,markerOption{:},'MarkerFaceColor', color_blue)  
        text(i,j,num2str(cnt),textOption{:})
        cnt=cnt+1;
    end
end
set(gca,axesOption{:})
axis xy % <- have bottom-left part be the origin
%% if 3rd argument (anomalous node indices) is given, plot them in red
if nargin==3
    % show red seeds
%     [ix,iy]= ind2sub([nx,ny],idx_anom);         
%     plot(ix,iy,markerOption{:},'MarkerFaceColor', color_red)    
    for i=1:length(idx_anom)
        [ix,iy]= ind2sub([nx,ny],idx_anom(i));             
        plot(ix,iy,markerOption{:},'MarkerFaceColor', color_red)    
        text(ix,iy,num2str(idx_anom(i)),textOption{:})
    end
    set(gca,axesOption{:})
end

drawnow
% keyboard