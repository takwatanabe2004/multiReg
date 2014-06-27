%% save_347_graphinfo_2d.m
% (05/28/2014)
%-------------------------------------------------------------------------%
% t_get_347_graphinfo_2d.m (01/07/2014)
% - Code from t_get_347_graphinfo.m....modified for 2-d slice for simulation
%--------------------------------------------------------------------------
% Code modified from t_create_6d_NN_adjmat.m...here we negate the sign of the 
% x-coordinate in the roiMNI...doing this will make the sampling order agree with
% the order of lexigraphic indexing (see t_j24_study_roiMNI_indexing2.m)
%--------------------------------------------------------------------------
%%
clear
purge

fsave=true;
% fsave=false;

load([get_rootdir,'/simulation/15by15/Results_tak.mat'], 'coord_used')
plot_option={'linewidth',3};
%%
outPath=[get_rootdir,'/simulation/15by15/graph_info347_2d.mat'];
outVars={'adjmat', 'coord', 'roiMNI', 'roiMNI_flipx', 'timeStamp', 'mFileName'};
%%
roiMNI=coord_used(:,1:2);

%==========================================================================
% The original index order of the roiMNI doesn't follow the convention of matlab....
%  - so flip the sign of the x-coordinate
%==========================================================================
roiMNI_flipx= bsxfun(@times, roiMNI, [-1 1]);
% return
%%
%==========================================================================
% - r: indexing of the 2d coordinate
% - rlex: lexicographic index order for the 2d coordinates r=(x,y)
%==========================================================================
[coord.r,coord.SPACING]=tak_discretize_coord_2d(roiMNI_flipx);
coord.num_nodes=size(coord_used,1);
coord.nx=max(coord.r(:,1));
coord.ny=max(coord.r(:,2));
coord.NSIZE=[coord.nx, coord.ny];
coord.N=prod(coord.NSIZE);
% return

% lexicographic ordering of r=(x,y)
coord.rlex = sub2ind(coord.NSIZE, coord.r(:,1), coord.r(:,2));

% check if the seeds are sampled in lexicographic order
isequal(coord.rlex,sort(coord.rlex))
% return
%% create 6d-coordinate...which are coordinates of the edges
%==========================================================================
% s: 4d coordinates... s=(r1,r2), r1=(x1,y1), r2=(x2,y2)
% l: index pair of 2-d coordinates, where the 2d coordinates are 
%    represented lexicographically
%==========================================================================
p = nchoosek(coord.num_nodes,2);
coord.s = zeros(p,4);
coord.l = zeros(p,2);
cnt=1;

for jj=1:coord.num_nodes
    for ii=jj+1:coord.num_nodes
        coord.s(cnt,:)=[coord.r(ii,:),coord.r(jj,:)];
        coord.l(cnt,1)=coord.rlex(ii);
        coord.l(cnt,2)=coord.rlex(jj);
        cnt=cnt+1;
    end
end

if ~all(coord.l(:,1)>coord.l(:,2))
    error('messed up again...')
end
% return
%==========================================================================
% lexicograph ordering of the 6d coordinates
%==========================================================================
coord.slex=sub2ind([coord.NSIZE,coord.NSIZE], ...
    coord.s(:,1),coord.s(:,2),...
    coord.s(:,3),coord.s(:,4));

figure,imexp
subplot(1,3,[1,2]),tplot(coord.r),legend('x','y','z')
subplot(133),tplot(coord.rlex)
figure,imexp
subplot(2,3,[1,2]),tplot(coord.s(:,1:2))
legend('x','y'),xlim([1,1000])
subplot(2,3,[4,5]),tplot(coord.s(:,3:4))
legend('xx','yy')
subplot(2,3,[3,6]),tplot(coord.slex)
tplott(coord.l)

% check if the sampling order agrees with the lexicograph order
isequal(coord.slex,sort(coord.slex))
% return
%% create a nearest-neighbor graph in 4d...do it brute force
adjmat=sparse(p,p);
timeTotal=tic;
tic
for i=1:p    
    % find the index-set of the 1st order nearest neighbor
    idx_NN=sum(abs(bsxfun(@minus, coord.s, coord.s(i,:))),2)==1;

    % nearest-neighbors
    adjmat(i,idx_NN)=1;
    
    %======================================================================
    % sanity check of the nearest neighbor property...
    % - the selected neighboring edge-sets should have l1 distance of 1 with the edge in question
    %======================================================================
    dist_set=abs(sum(bsxfun(@minus,coord.s(idx_NN==1,:),coord.s(i,:)),2));
    if ~all(bsxfun(@eq, dist_set, 1))
        error('argh..')
    end
end
timeTotal=toc(timeTotal);
timeStamp=tak_timestamp
mFileName=mfilename;
coord=orderfields(coord);

if fsave
    save(outPath,outVars{:})
end