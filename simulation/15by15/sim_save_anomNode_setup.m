% sim_save_anomNode_setup.m (01/07/2014)
% - figure out anomalous node configuration and save on disk
%==========================================================================
%--------------------------------------------------------------------------
%%
clear all
purge
load([get_rootdir,'/simulation/15by15/graph_info347_2d.mat'], 'coord')

fsave=true;
% fsave=false;

p=nchoosek(coord.num_nodes,2);
%%
outPath=[get_rootdir,'/simulation/15by15/sim_anom_node_info.mat'];
outVars={'anom_nodes','timeStamp','mFileName'};
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
%% figure out where to place the anomalous node clusters
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

%==========================================================================
% introduce a pair of clover-shaped anomalous node cluseters
%==========================================================================
% centroid of the anomalous nodes
coord_anom_ctr1=[4,3];
coord_anom_ctr2=[6,7];

%==========================================================================
% get index info for the 1st cluster of anom nodes
%==========================================================================
% index list of the anomlaous node locations
tmp1=[coord_anom_ctr1(1),coord_anom_ctr1(2)];    
%--- 4 nearest neighbors ---%
tmp1a=[coord_anom_ctr1(1)-1,coord_anom_ctr1(2)];
tmp1b=[coord_anom_ctr1(1)+1,coord_anom_ctr1(2)];
tmp1c=[coord_anom_ctr1(1),coord_anom_ctr1(2)-1];
tmp1d=[coord_anom_ctr1(1),coord_anom_ctr1(2)+1];

% (x,y) coordintae list of the anomlaous node locations
coord_anom1=[tmp1;tmp1a;tmp1b;tmp1c;tmp1d];
    
%--- get indices of the anomaluos nodes ---%
idx1=find(sum(bsxfun(@eq, coord.r, tmp1),2)==2);
idx1a=find(sum(bsxfun(@eq, coord.r, tmp1a),2)==2);
idx1b=find(sum(bsxfun(@eq, coord.r, tmp1b),2)==2);
idx1c=find(sum(bsxfun(@eq, coord.r, tmp1c),2)==2);
idx1d=find(sum(bsxfun(@eq, coord.r, tmp1d),2)==2);
coord_anom_ind1=[idx1;idx1a;idx1b;idx1c;idx1d];

%--- plot anomalous nodes in different color ---%
plot(tmp1(1),tmp1(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp1a(1),tmp1a(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp1b(1),tmp1b(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp1c(1),tmp1c(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp1d(1),tmp1d(2),markerOption{:},'MarkerFaceColor', color_red)
text(tmp1(1),tmp1(2),num2str(idx1),textOption{:})
text(tmp1a(1),tmp1a(2),num2str(idx1a),textOption{:})
text(tmp1b(1),tmp1b(2),num2str(idx1b),textOption{:})
text(tmp1c(1),tmp1c(2),num2str(idx1c),textOption{:})
text(tmp1d(1),tmp1d(2),num2str(idx1d),textOption{:})
coord_anom1
% coord_anom_ind1=sort(coord_anom_ind1)
% pause
%==========================================================================
% get index info for the 2nd cluster of anom nodes
%==========================================================================
% index list of the anomlaous node locations
tmp2=[coord_anom_ctr2(1),coord_anom_ctr2(2)];    
%--- 4 nearest neighbors ---%
tmp2a=[coord_anom_ctr2(1)-1,coord_anom_ctr2(2)];
tmp2b=[coord_anom_ctr2(1)+1,coord_anom_ctr2(2)];
tmp2c=[coord_anom_ctr2(1),coord_anom_ctr2(2)-1];
tmp2d=[coord_anom_ctr2(1),coord_anom_ctr2(2)+1];

% (x,y) coordintae list of the anomlaous node locations
coord_anom2=[tmp2;tmp2a;tmp2b;tmp2c;tmp2d];
    
%--- get indices of the anomaluos nodes ---%
idx1=find(sum(bsxfun(@eq, coord.r, tmp2),2)==2);
idx1a=find(sum(bsxfun(@eq, coord.r, tmp2a),2)==2);
idx1b=find(sum(bsxfun(@eq, coord.r, tmp2b),2)==2);
idx1c=find(sum(bsxfun(@eq, coord.r, tmp2c),2)==2);
idx1d=find(sum(bsxfun(@eq, coord.r, tmp2d),2)==2);
coord_anom_ind2=[idx1;idx1a;idx1b;idx1c;idx1d];

%--- plot anomalous nodes in different color ---%
plot(tmp2(1),tmp2(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp2a(1),tmp2a(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp2b(1),tmp2b(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp2c(1),tmp2c(2),markerOption{:},'MarkerFaceColor', color_red)
plot(tmp2d(1),tmp2d(2),markerOption{:},'MarkerFaceColor', color_red)
text(tmp2(1),tmp2(2),num2str(idx1),textOption{:})
text(tmp2a(1),tmp2a(2),num2str(idx1a),textOption{:})
text(tmp2b(1),tmp2b(2),num2str(idx1b),textOption{:})
text(tmp2c(1),tmp2c(2),num2str(idx1c),textOption{:})
text(tmp2d(1),tmp2d(2),num2str(idx1d),textOption{:})
    %     return
% coord_anom_ind2=sort(coord_anom_ind2)
%%
%==========================================================================
% figure out the index set corresponding to the connections among the
% anomalous node cluster pairs
%==========================================================================
% do it by brute force loop (for sanity check)
anom_conn_idx=[];
cnt=1;
for i=1:coord.num_nodes
    flag1=sum(bsxfun(@eq, coord_anom_ind1, i))==1;
%     if sum(bsxfun(@eq, coord_anom_ind1, i))==1
%         i
%     end
    for j=i+1:coord.num_nodes
%         if sum(bsxfun(@eq, coord_anom_ind2, j))==1
%             j
%         end
        flag2=sum(bsxfun(@eq, coord_anom_ind2, j))==1;
        
        if flag1 && flag2
            anom_conn_idx=[anom_conn_idx;cnt];
        end
        cnt=cnt+1;
    end
end

% a more elegant way
TMP = tak_dvecinv(1:p,0);
anom_conn_idx2=[];
for i=1:length(coord_anom_ind1)
    anom1=coord_anom_ind1(i);
    for j=1:length(coord_anom_ind2)        
        anom2=coord_anom_ind2(j);
        anom_conn_idx2=[anom_conn_idx2;TMP(anom1,anom2)];
    end
end

if isequal(anom_conn_idx,sort(anom_conn_idx2))
    disp('----- good! sanity check complete -----')
else
    error('----- welp...more debugging... -----')    
end
%%
%==========================================================================
% mask representing the locations of the anomalous nodes in matrix from
%==========================================================================
anom_mask = false(coord.num_nodes);
for i=1:length(coord_anom_ind1)
    ii=coord_anom_ind1(i);
    for j=1:length(coord_anom_ind2)
        jj=coord_anom_ind2(j);
        anom_mask(ii,jj)=true;
    end
end
% symmetrize
anom_mask = anom_mask + anom_mask';
anom_mask_vec = tak_dvec(anom_mask);
imedgel(anom_mask)
tplott(anom_mask_vec)
% return
%%

% pool relevant info into a struct
anom_nodes.coor_ctr1=coord_anom_ctr1;
anom_nodes.coord1=coord_anom1;
anom_nodes.coordind1=coord_anom_ind1;

anom_nodes.coor_ctr2=coord_anom_ctr2;
anom_nodes.coord2=coord_anom2;
anom_nodes.coordind2=coord_anom_ind2;

anom_nodes.idx_conn = anom_conn_idx
anom_nodes.mask = anom_mask;
anom_nodse.maskVec=anom_mask_vec;

anom_nodes.coordind1
anom_nodes.coordind2

if fsave
    timeStamp=tak_timestamp;
    mFileName=mfilename;
    save(outPath,outVars{:})
end