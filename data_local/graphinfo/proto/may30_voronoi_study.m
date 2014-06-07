%% may30_voronoi_study.m
% (05/30/2014)
%=========================================================================%
% Study voronoi diagram properties....this could be useful for defining
% the nearest-neighbor graph for the WashU parcellation
%=========================================================================%
clear
purge

x = gallery('uniformdata',[1 10],0);
y = gallery('uniformdata',[1 10],1);
voronoi(x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% food for thought: define nearest-neighbor based on Voronoi set?
%   - differences among adjacent cells well be penalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%