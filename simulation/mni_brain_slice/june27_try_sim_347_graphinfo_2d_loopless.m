%% june27_try_sim_347_graphinfo_2d_loopless.m
% (06/27/2014)
%=========================================================================%
% - Try to replicate sim_save_347_graphinfo_2d.m without resorting to loops
%=========================================================================%
%%
clear all;
purge

rootdir = fileparts(mfilename('fullpath'));
load([rootdir,'/Results_tak.mat'], 'coord_used')
plot_option={'linewidth',3};
%% cell blocks from sim_save_347_graphinfo_2d.m
roiMNI=coord_used(:,1:2);

%==========================================================================
% The original index order of the roiMNI doesn't follow the convention of matlab....
%  - so flip the sign of the x-coordinate
%==========================================================================
roiMNI_flipx= bsxfun(@times, roiMNI, [-1 1]);
% return

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

% lexicographic ordering of r=(x,y)
coord.rlex = sub2ind(coord.NSIZE, coord.r(:,1), coord.r(:,2));

% check if the seeds are sampled in lexicographic order
if ~isequal(coord.rlex,sort(coord.rlex))
    error('messed up again...')
end


%%%%% create 6d-coordinate...which are coordinates of the edges %%%%%%%%%%%
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

%==========================================================================
% lexicograph ordering of the 6d coordinates
%==========================================================================
coord.slex=sub2ind([coord.NSIZE,coord.NSIZE], ...
    coord.s(:,1),coord.s(:,2),...
    coord.s(:,3),coord.s(:,4));

% figure,imexp
% subplot(1,3,[1,2]),tplot(coord.r),legend('x','y','z')
% subplot(133),tplot(coord.rlex)
% figure,imexp
% subplot(2,3,[1,2]),tplot(coord.s(:,1:2))
% legend('x','y'),xlim([1,1000])
% subplot(2,3,[4,5]),tplot(coord.s(:,3:4))
% legend('xx','yy')
% subplot(2,3,[3,6]),tplot(coord.slex)
% tplott(coord.l)

% check if the sampling order agrees with the lexicograph order
if ~isequal(coord.slex,sort(coord.slex))
    error('messed up again...')
end
%% from sim_save_347_graphinfo_2d.m: create a NN-graph in 4d by brute force
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


timeTotal=toc(timeTotal)
%% loopless way? sample edge-entries from full diffmat
nx=coord.nx;
ny=coord.ny;
d=nx*ny; % # <- # nodes

C4d=tak_diffmat_4d( [coord.nx, coord.ny, coord.nx, coord.ny],0);

%-------------------------------------------------------------------------%
% any interesting thing going on with inc2adj?
%-------------------------------------------------------------------------%
adj4d = inc2adj(C4d);
Lap4d=C4d'*C4d;
figure,imexpb
subplot(131),tspy(C4d)
subplot(132),tspy(Lap4d)
subplot(133),tspy(adj4d)

% d^2*4 % <- # edges on C4d_circ
% size(C4d,1)/size(C4d_brute,1)

%-------------------------------------------------------------------------%
% create a d^2-length mask-vector indicating 4-d coordinates included
% in the irregular 347-node coordinate
%-------------------------------------------------------------------------%
mask_vec=false(d^2,1);
mask_vec(coord.slex)=true;
% figure,imexpb,tplot(mask_vec)

%-------------------------------------------------------------------------%
% sample adjmat indices cooresponding to the above mask
%-------------------------------------------------------------------------%
adj4d_subsamp = adj4d(coord.slex,coord.slex);
adj4d_subsamp = adj4d_subsamp + adj4d_subsamp'; % symmetrie
C4d_subsamp = tak_adjmat2incmat(adj4d_subsamp);
Lap4d_subsamp=C4d_subsamp'*C4d_subsamp;
figure,imexpb
subplot(131),tspy(C4d_subsamp)
subplot(132),tspy(Lap4d_subsamp)
subplot(133),tspy(adj4d_subsamp)

%=========================================================================%
% compare with the brute-force method
%=========================================================================%
C4d_brute = tak_adjmat2incmat(adjmat);
Lap4d_brute=C4d_brute'*C4d_brute;
figure,imexpb
subplot(131),tspy(C4d_brute)
subplot(132),tspy(Lap4d_brute)
subplot(133),tspy(adjmat)

isequal(C4d_brute,C4d_subsamp)
isequal(adjmat,adj4d_subsamp)
% tspyl(adjmat)
% tspyl(inc2adj(C4d_brute))
% tspyl(inc2adj(C4d_brute)+inc2adj(C4d_brute)')