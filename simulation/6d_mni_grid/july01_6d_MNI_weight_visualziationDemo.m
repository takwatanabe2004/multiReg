%% july01_6d_MNI_weight_visualziationDemo
% (07/01/2014)
%=========================================================================%
% - Demo for some "convenience" scripts i wrote to visualize weight vector results
%=========================================================================%
%%
clear
purge

GRID='Grid326'; % {'Grid326','Grid1068','WashU'}

rootdir = fileparts(mfilename('fullpath'));

randn('state',0)
rand('state',0)
%% load ground-truth weight support and graphInfo of node parcellation
dataPath = [rootdir,'/anomNodeCluster2','_',GRID,'.mat'];
dataVars = {'wsupp', 'idx_anomNodeClusters'};
load(dataPath,dataVars{:})
idx_supp = find(wsupp);

% load graphInfo
dataPath = [get_rootdir,'/data_local/graphinfo/graph_info_',GRID,'.mat'];
dataVars = {'C','coord'};
load(dataPath,dataVars{:})

% feature info
d = prod(coord.nsamp);
p = nchoosek(d,2);

% circmat for fft-based FL
NSIZE = [coord.NSIZE, coord.NSIZE];
Ccirc = tak_diffmat(NSIZE,1);
[A,b] = tak_get_augMat_maskMat(NSIZE,coord.slex);
%% visualization demo (suppress later)
% figure,imexpb
% subplot(121),tspy(C),axis normal
% subplot(122),tspy(C'*C)
% figure,imexpb
% subplot(121),tspy(Ccirc),axis normal
% subplot(122),tspy(Ccirc'*Ccirc)
% drawnow


imconnEdgel(wsupp)
imconnEdger(wsupp)
disp('hit key'); pause

%=========================================================================%
w=zeros(size(wsupp));
idx=find(wsupp);
w(idx)=rand(size(idx))-.5;

%=========================================================================%
purge
figure,imconnYeo(w,GRID,1)
figure,imconnYeo(w,GRID,1,1)
imconnYeol(w,GRID,1)
imconnYeol(w,GRID,1,1)
imconnYeor(w,GRID,1,1)
imconnYeor(w,GRID,1)
disp('hit key'); pause

%=========================================================================%
purge
figure,imconnYeoEdge(w,GRID)
figure,imconnYeoEdge(w,GRID,1)
imconnYeoEdgel(w,GRID)
imconnYeoEdgel(w,GRID,1)
imconnYeoEdger(w,GRID)
imconnYeoEdger(w,GRID,1)
disp('hit key'); pause

%=========================================================================%
purge
figure,imconnYeoTriag(w,GRID)
figure,imconnYeoTriag(w,GRID,1)
imconnYeoTriagl(w,GRID)
imconnYeoTriagl(w,GRID,1)
imconnYeoTriagr(w,GRID)
imconnYeoTriagr(w,GRID,1)