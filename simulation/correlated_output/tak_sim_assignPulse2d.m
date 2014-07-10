function w = tak_sim_assignPulse2d(NSIZE, nPatches)
% w = tak_sim_assignPatch1d(p, nPatches)
%=========================================================================%
% - Assign "nPatches" of "nonoverlapping" patches on a 1d signal
% (patchLenBase = "base-length" of each patches)
%=========================================================================%
% (07/10/2014)
%%
nx=NSIZE(1);
ny=NSIZE(2);
W=zeros(nx,ny);

radList = round(min(nx,ny)./[5:8]);

[xx,yy]=ndgrid(1:nx,1:ny);

for ipatch = 1:nPatches
    %=====================================================================%
    % center points and radius
    %=====================================================================%
    rad = radList(randsample(length(radList),1));
    
    %---------------------------------------------------------------------% 
    % select center point in such a way that elliptical signal won't spill 
    % out the FOV
    %---------------------------------------------------------------------% 
    ctrx = randsample(nx,1);
    ctry = randsample(ny,1);
    ctr = [ctrx,ctry];
    
    %=====================================================================%
    % get (x,y) coordinates of the cluster
    %=====================================================================% 
    quad_kern = @(crd) (1-crd.^2).^2 .* (abs(crd)<=1);
    crd = @(ctr,rad) sqrt( ((xx-ctr(1))/rad).^2 + ((yy-ctr(2))/rad).^2  + ...
                            2*(rand-.5)*(xx-ctr(1)).*(yy-ctr(2))/rad^2  );

    %=====================================================================%
    % arbitrary magnitude for the cluster
    %=====================================================================%
    mag =  tak_sample_signed_unif([3,5],1);

    %=====================================================================%
    % final support of the weigth vector in 2d space
    %=====================================================================%    
    tmp = zeros(nx,ny);
%     keyboard
    tmp = tmp + mag*quad_kern(crd(ctr,rad));
%     keyboard
    W = W + tmp;
end
w=W(:);