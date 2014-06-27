%% plot_node_orientation.m
% (06/26/2014)
%==========================================================================
% (08/29/2013) 
% - Create figures for showing the seed placement for the simulation, with the
%   'healthy' representing all blue nodes, and 'diseaes' represented with few
%   red nodes
% - This script has to be ran before sim_samples.m (see last cell in this code)
%===========================================================================
clear
purge

rootdir=fileparts(mfilename('fullpath'));
%% setup output figname
figname1=[rootdir,'/sim_seeds_HC'];
figname2=[rootdir,'/sim_seeds_DS'];

% fsave=true;
fsave=false;
%% load data
dataPath1 = [rootdir,'/../sim_data_dataset_setup7.mat']
dataVars1={'sim_opt','DISTR'};
load(dataPath1,dataVars1{:})

NSIZE=sim_opt.NSIZE;
ix_anom=sim_opt.ix_anom;
iy_anom=sim_opt.iy_anom;
%%
screenSize = [10 100 900 900];
msize=28;
mwidth=4;

% set background color
colorBack=[0.85 1 0.85];

axesOption={'XTick',[1,5,10,15],'YTick',[1,5,10,15],'linewidth',3,'Fontweight','b',...
    'fontsize',36,'TickLength',[0 0],'xaxislocation','bottom'};

% color_blue = [0.8 0.8 1];
% color_red = [1, 0.8 0.8];
color_blue = [0.6 0.6 1];
color_red = [1, 0.6 0.6];
markerOption={'o', 'MarkerEdgeColor','k','MarkerSize',msize,'linewidth',mwidth};
%%
% superimposing 2d plots over an image is the easiest way to go for me
% figure('visible','off')
figure,set(gcf,'Units','pixels','Position', screenSize)
imagesc(zeros(NSIZE))

% set background color
colormap(colorBack) 

% show blue seeds first
hold on
axis('on','image');
for i=1:NSIZE(1)
    for j=1:NSIZE(2)
        plot(i,j,markerOption{:},'MarkerFaceColor', color_blue)  
    end
end
set(gca,axesOption{:})

% save as healthy seeds
if fsave
    savefig(figname1,'pdf')
%     savefig(figname1,'png')
end

% show red seeds
for i=1:length(ix_anom)
    for j=1:length(iy_anom)
        ix = ix_anom(i);
        iy = iy_anom(j);
        plot(ix,iy,markerOption{:},'MarkerFaceColor', color_red)    
        plot(iy,ix,markerOption{:},'MarkerFaceColor', color_red)    
    end
end
set(gca,axesOption{:})

if fsave
    savefig(figname2,'pdf')
%     savefig(figname2,'png')
end