%% june11_make_fake_target.m
% (06/11/2014)
%=========================================================================%
% Make artificial target data (y) for regression
%=========================================================================%
%%
clear
purge

grid='Grid326'; % {'Grid326','Grid1068','WashU'}

%% configs
snr = 0;
seedPoints.randn=0;
seedPoints.rand=0;
rand('state',seedPoints.rand)
randn('state',seedPoints.randn)

outPath = [fileparts(mfilename('fullpath')),...
    '/fake_target_',grid,'_snr',num2str(snr),'.mat']
outVars={'y','wtrue','snr','seedPoints','mFileName','timeStamp'};
% return
%%
load([get_rootdir,'/data_local/designMatrix_FC_',grid,'.mat'],'X')

%-------------------------------------------------------------------------%
% get yeo info
%-------------------------------------------------------------------------%
load([get_rootdir,'/data_local/yeoLabelInfo/yeo_info_',grid,'_dilated5mm.mat'],...
    'roiLabel','yeoLabels')

[n,p]=size(X);
d=length(roiLabel);
%% create mask indicating the support of the artificial weight vector
%=========================================================================%
% idxsrt:       groups the nodes according to yeoPlus network
% idxsrt_back:  reverse mapping of idxsrt
%=========================================================================%
[idxsrt,labelCount] = tak_get_yeo_sort(roiLabel);
[~,idxsrt_back] = sort(idxsrt);
% idxsrt(idxsrt_back)
% return

%=========================================================================%%
% the yeo list:
%-------------------------------------------------------------------------%
%  1   'Visual'
%  2   'Somatomotor'
%  3   'Dorsal Attention'
%  4   'Ventral Attention'
%  5   'Limbic'
%  6   'Frontoparietal'
%  7   'Default'
%  8   'Striatum'
%  9   'Amygdala'
% 10   'Hippocampus'
% 11    'Thalamus'
% 12    'Cerebellum'
% 13    'Unlabeled'
%=========================================================================%
% impacted networks 
% netList = [6 1; 4 3; 7 7; 7 4; 5 5];
netList = [6 1; ];

nLabel = length(yeoLabels);
roiLabel_srt = roiLabel(idxsrt);

mask_mat = false(d);
for i=1:size(netList,1)
    mask_mat(roiLabel_srt==netList(i,1),roiLabel_srt==(netList(i,2)))=true;
end
mask_mat = tril(mask_mat,-1); % <- just need the lower triangular part
% imedgel(mask_mat),  axis off,    %colorbar off%colorbar('location','northoutside')
%     tak_local_linegroups6_nocirc(gcf,labelCount,yeoLabels)
% return
%% create "fake" weight vector by adding values on the support/mask 
mask_vec=tak_dvec(mask_mat);
wtrue_srt = zeros(p,1);
wtrue_srt(mask_vec)=2*rand(sum(mask_vec),1) -1; % <- uniform from [-1,1]
wtrue_srt_mat = tak_dvecinv(wtrue_srt,0);
    
%-------------------------------------------------------------------------%
% map weight vector from "sorted" space to its native space
%-------------------------------------------------------------------------%
wtrue_mat = wtrue_srt_mat(idxsrt_back,idxsrt_back);
wtrue = tak_dvec(wtrue_mat);
% imcovvl(wtrue_srt_mat),  axis off,    colorbar('location','northoutside')
%     tak_local_linegroups6_nocirc(gcf,labelCount,yeoLabels)
% imcovvl(wtrue_mat)
% imcovvl(wtrue_mat(idxsrt,idxsrt))
% tplott(wtrue),imexpt
% return
%% make fake data via additive noise model: y=X*w + noise;
snr=0;
X = randn(n,p);
y = X*wtrue + snr*randn(n,1);

%% save
mFileName = mfilename('fullpath');
timeStamp=tak_timestamp;
save(outPath,outVars{:})