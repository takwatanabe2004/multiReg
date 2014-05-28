function [rdiscr, dr] = tak_discretize_coord_2d(r)
%| rd = tak_discretize_coord_2d(r)
%|--------------------------------------------------------------------------------|
%|      r=(x,y)-coordinates in real domain
%| rdiscr=(xd,yd)- discretized coordinates {1,...,nx},{1,...,ny}
%|     dr=(dx,dy)-spacing in the x,y direction
%|
%| ...method sorta specific for the slab data (not for general 2d data)
%|--------------------------------------------------------------------------------|
%|        r:  (2 x ncoord)
%|   rdiscr: (2 x ncoord)
%|   dr
%|--------------------------------------------------------------------------------|
%| (05/21/2013)
%%
x=r(:,1);
y=r(:,2);

xx=unique(x);
yy=unique(y);

dx=xx(2)-xx(1);
dy=yy(2)-yy(1);

x=(x-min(x))/dx+1;
y=(y-min(y))/dy+1;
rdiscr=[x,y];
dr=[dx,dy];