
clear
purge

load([get_rootdir,'/simulation/15by15/graph_info347_2d.mat'], 'coord')
% load('sim_anom_node_info.mat', 'anom_nodes')
load('sim_anom_node_info_8nn.mat', 'anom_nodes')
%% set figure options
msize=28;
mwidth=4;

% set background color
colorBack=[0.85 1 0.85];

% marker colors
color_blue = [0.6 0.6 1];
color_red = [1, 0.6 0.6];

axesOption={'Xtick',[1:coord.NSIZE(1)],'Ytick',[1:coord.NSIZE(2)],...
    'TickLength',[0 0],'xaxislocation','bottom','fontsize',14,'fontweight','b','box','on'};
markerOption={'o', 'MarkerEdgeColor','k','MarkerSize',msize,'linewidth',mwidth};
textOption={'fontsize',16,'fontweight','b','HorizontalAlignment','center',...
            'VerticalAlignment','middle'};
%%
figure,set(gcf,'Units','pixels','Position', [210 100 600 600]),hold on
imagesc(zeros(coord.NSIZE)')
xlabel('x',textOption{:})
ylabel('y',textOption{:})
% set background color
colormap(colorBack) 
for i=1:coord.num_nodes
    plot(coord.r(i,1),coord.r(i,2),markerOption{:},'MarkerFaceColor', color_blue)  
    text(coord.r(i,1),coord.r(i,2),num2str(i),textOption{:})
end
axis image
set(gca,axesOption{:})
drawnow

%% plot anom nodes
% purge
%---  first plot healthy nodes ---%
figure,set(gcf,'Units','pixels','Position', [810 100 600 600]),hold on
imagesc(zeros(coord.NSIZE)')
xlabel('x',textOption{:})
ylabel('y',textOption{:})
% set background color
colormap(colorBack) 
for i=1:coord.num_nodes
    plot(coord.r(i,1),coord.r(i,2),markerOption{:},'MarkerFaceColor', color_blue)  
    text(coord.r(i,1),coord.r(i,2),num2str(i),textOption{:})
end
axis image
set(gca,axesOption{:})
drawnow

%--- plot anomalous nodes in different color ---%
% first cluster
for i=1:length(anom_nodes.coord1)
    x=anom_nodes.coord1(i,1);
    y=anom_nodes.coord1(i,2);
    plot(x,y,markerOption{:},'MarkerFaceColor', color_red)
    text(x,y,num2str(anom_nodes.coordind1(i)),textOption{:})
end

% second cluster
for i=1:length(anom_nodes.coord2)
    x=anom_nodes.coord2(i,1);
    y=anom_nodes.coord2(i,2);
    plot(x,y,markerOption{:},'MarkerFaceColor', color_red)
    text(x,y,num2str(anom_nodes.coordind2(i)),textOption{:})
end

% plot(tmp1(1),tmp1(2),markerOption{:},'MarkerFaceColor', color_red)
% plot(tmp1a(1),tmp1a(2),markerOption{:},'MarkerFaceColor', color_red)
% plot(tmp1b(1),tmp1b(2),markerOption{:},'MarkerFaceColor', color_red)
% plot(tmp1c(1),tmp1c(2),markerOption{:},'MarkerFaceColor', color_red)
% plot(tmp1d(1),tmp1d(2),markerOption{:},'MarkerFaceColor', color_red)
% text(tmp1(1),tmp1(2),num2str(idx1),textOption{:})
% text(tmp1a(1),tmp1a(2),num2str(idx1a),textOption{:})
% text(tmp1b(1),tmp1b(2),num2str(idx1b),textOption{:})
% text(tmp1c(1),tmp1c(2),num2str(idx1c),textOption{:})
% text(tmp1d(1),tmp1d(2),num2str(idx1d),textOption{:})
% coord_anom1
% coord_anom_ind1=sort(coord_anom_ind1)