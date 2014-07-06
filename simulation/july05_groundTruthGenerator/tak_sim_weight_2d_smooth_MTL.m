function [W,idx_supp] = tak_sim_weight_2d_smooth_MTL(nx,ny,q)
% [W,idx_supp] = tak_sim_weight_2d_smooth_MTL(nx,ny,q)
%=========================================================================%
% - Generate 2d grouth truth weight vector.
%=========================================================================%
% (07/05/2014)
%%
p=nx*ny;
[xx,yy]=ndgrid(1:nx,1:ny);

%=========================================================================%
% create 3 clusters specified "center points" and "radius"
%=========================================================================%
ctr1 = round([nx*.5 ny*.2]);
rad1 = round( min(nx,ny)/10 );

ctr2 = round([nx*.3 ny*.6]);
rad2 = round( min(nx,ny)/6 );

ctr3 = round([nx*.75 ny*.75]);
rad3 = round( min(nx,ny)/8 );
% ctr1 = [21 61]; rad1 = 10;
% ctr2 = [66 66]; rad2 = 8;
% ctr3 = [44 25]; rad3 = 18;

%=========================================================================%
% get (x,y) coordinates of the 3 clusters
%=========================================================================% 
quad_kern = @(crd) (1-crd.^2).^2 .* (abs(crd)<=1);
crd1 = @(ctr,rad) sqrt( (xx-ctr(1)).^2 + (yy-ctr(2)).^2 )/rad;
% crd = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2 );
crd2 = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2  + ...
                        1*(xx-ctr(1)).*(yy-ctr(2))/rad^2               );
crd3 = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2  + ...
                        -1*(xx-ctr(1)).*(yy-ctr(2))/rad^2               );
                    
W = zeros(p,q);
for i=1:q
    % random magnitude for each cluster
    mag =  tak_sample_signed_unif([3,5],3);

    %=====================================================================%
    % final support of the weigth vector in 2d space
    %=====================================================================%
    Wtmp = zeros(nx,nx);
    Wtmp = Wtmp + mag(1)*quad_kern(crd1(ctr1,rad1)) ...
                + mag(2)*quad_kern(crd2(ctr2,rad2)) ...
                + mag(3)*quad_kern(crd3(ctr3,rad3));
            
    W(:,i) = Wtmp(:);
end
idx_supp = find(W(:,1));
