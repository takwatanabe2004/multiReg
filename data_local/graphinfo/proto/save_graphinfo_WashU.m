%% save_graphinfo_WashU.m
% (05/30/2014)
%=========================================================================%
% - get (p x p) adjacency matrix [adjmat] and incidence matrix for the
%   WashU parcelation....which does NOT follow a "grid-based" parecellation
%--------------------------------------------------------------------------
% Code modified from t_create_6d_NN_adjmat.m...here we negate the sign of the 
% x-coordinate in the roiMNI...doing this will make the sampling order agree with
% the order of lexigraphic indexing (see t_j24_study_roiMNI_indexing2.m)
%=========================================================================%
%%
clear
purge

fsave=true;
parcellation = 'WashU'; 
%% 
outPath=[get_rootdir, '/data_local/graphinfo/graph_info_',parcellation,'.mat']
outVars={'adjmat', 'C', 'coord', 'roiMNI', 'roiMNI_flipx', 'timeStamp',...
    'timeTotal', 'mFileName'};
% return
%%
yeoInfoPath=[get_rootdir, ...
    '/data_local/yeoLabelInfo/yeo_info_',parcellation,'.mat'];
load(yeoInfoPath,'yeoLabels', 'roiMNI', 'roiLabel')
plot_option={'linewidth',3};

d=size(roiMNI,1);
p=nchoosek(d,2);
%% evaluate nearest-neighbor nodes 3-D space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - the "grid-based parcellations" nodes were lexigraphically ordered from
%   the beginning (after flipping the x-ccoordinate - which is why I 
%   defined "roiMNI_flipx" in save_graphinfo_grid.m)
% 
% - plus, WashU parcellation are in arbitrary coordinates in the MNI space, 
%   so more effort is need to evaluate the nearest neighborhood structure
%
% - the neighborhood structure will follow a "nearest-centroid" approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%=========================================================================%
% The original index order of the roiMNI doesn't follow matlab convention....
%  - so flip the sign of the x-coordinate
%=========================================================================%
roiMNI_flipx= bsxfun(@times, roiMNI, [-1 1 1]);

%=========================================================================%
% - r: indexing of the 3d coordinate
% - rlex: lexicographic index order for the 3d coordinates r=(x,y,z)
%=========================================================================%
[coord.r,coord.SPACING]=tak_discretize_coord(roiMNI_flipx);
coord.nsamp=d;
coord.nx=max(coord.r(:,1));
coord.ny=max(coord.r(:,2));
coord.nz=max(coord.r(:,3));
coord.NSIZE=[coord.nx, coord.ny, coord.nz];
coord.N=prod(coord.NSIZE);

% lexicographic ordering of r=(x,y,z)
coord.rlex = sub2ind(coord.NSIZE, coord.r(:,1), coord.r(:,2), coord.r(:,3));

% check if the seeds are sampled in lexicographic order
if ~isequal(coord.rlex,sort(coord.rlex))
    error('meh....debug')
end
%% create 6d-coordinate...which are coordinates of the edges
%==================================================================================
% s: 6d coordinates... s=(r1,r2), r1=(x1,y1,z1), r2=(x2,y2,z2)
% l: index pair of 3-d coordinates, where the 3d coordinates are represented lexicographically
%==================================================================================
coord.s = zeros(p,6);
coord.l = zeros(p,2);
cnt=1;

for jj=1:d
    for ii=jj+1:d
        coord.s(cnt,:)=[coord.r(ii,:),coord.r(jj,:)];
        coord.l(cnt,1)=coord.rlex(ii);
        coord.l(cnt,2)=coord.rlex(jj);
        cnt=cnt+1;
    end
end

if ~all(coord.l(:,1)>coord.l(:,2))
    error('messed up again...')
end

%==================================================================================
% lexicograph ordering of the 6d coordinates
%==================================================================================
coord.slex=sub2ind([coord.NSIZE,coord.NSIZE], ...
    coord.s(:,1),coord.s(:,2),coord.s(:,3), ...
    coord.s(:,4),coord.s(:,5),coord.s(:,6));

% figure,imexp
% subplot(1,3,[1,2]),tplot(coord.r),legend('x','y','z')
% subplot(133),tplot(coord.rlex)
% figure,imexp
% subplot(2,3,[1,2]),tplot(coord.s(:,1:3))
% legend('x','y','z'),xlim([1,1000])
% subplot(2,3,[4,5]),tplot(coord.s(:,4:6))
% legend('xx','yy','zz')
% subplot(2,3,[3,6]),tplot(coord.slex)
% 
% tplott(coord.l)


% check if the sampling order agrees with the lexicograph order
isequal(coord.slex,sort(coord.slex))
%% create a nearest-neighbor graph in 6d...do it brute force
adjmat=sparse(p,p);
timeTotal=tic;
tic
for i=1:p
    if mod(i,2000)==0; i,toc, end;
    
    % find the index-set of the 1st order nearest neighbor
    idx_NN=sum(abs(bsxfun(@minus, coord.s, coord.s(i,:))),2)==1;

    % nearest-neighbors
    adjmat(i,idx_NN)=1;
    
    %==================================================================================
    % sanity check of the nearest neighbor property...
    % - the selected neighboring edge-sets should have l1 distance of 1 with the edge in question
    %==================================================================================
    dist_set=abs(sum(bsxfun(@minus,coord.s(idx_NN==1,:),coord.s(i,:)),2));
    if ~all(bsxfun(@eq, dist_set, 1))
        error('argh..')
    end
end

%-------------------------------------------------------------------------%
% incidence matrix (aka finite differencing matrix
%-------------------------------------------------------------------------%
C = tak_adjmat2incmat(adjmat);

timeTotal=toc(timeTotal);
timeStamp=tak_timestamp
mFileName=mfilename('fullpath');
coord=orderfields(coord);

if fsave
    save(outPath,outVars{:})
end
% toc