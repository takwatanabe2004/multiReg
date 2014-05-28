function [rdiscr, dr] = tak_discretize_coord(r)
%| rd = tak_discretize_coord(r)
%|--------------------------------------------------------------------------------|
%|      r=(x,y,z)-coordinates in real domain
%| rdiscr=(xd,yd,zd)- discretized coordinates {1,...,nx},{1,...,ny},{1,...,nz}
%|     dr=(dx,dy,dz)-spacing in the x,y,z direction
%|
%| ...method sorta specific for the slab data (not for general 3d data)
%|--------------------------------------------------------------------------------|
%|        r:  (3 x ncoord)
%|   rdiscr: (3 x ncoord)
%|   dr
%|--------------------------------------------------------------------------------|
%| (05/21/2013)
%%
x=r(:,1);
y=r(:,2);
z=r(:,3);

xx=unique(x);
yy=unique(y);
zz=unique(z);

dx=xx(2)-xx(1);
dy=yy(2)-yy(1);
dz=zz(2)-zz(1);

x=(x-min(x))/dx+1;
y=(y-min(y))/dy+1;
z=(z-min(z))/dz+1;
rdiscr=[x,y,z];
dr=[dx,dy,dz];