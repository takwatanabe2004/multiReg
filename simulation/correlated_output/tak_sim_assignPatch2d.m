function w = tak_sim_assignPatch2d(NSIZE, nPatches,patchLenBase)
% w = tak_sim_assignPatch1d(p, nPatches)
%=========================================================================%
% - Assign "nPatches" of "nonoverlapping" patches on a 1d signal
% (patchLenBase = "base-length" of each patches)
%=========================================================================%
% (07/10/2014)
%%
nx=NSIZE(1);
ny=NSIZE(2);
p = prod(NSIZE);
W=zeros(nx,ny);

if nargin==2
    patchLenBaseX= round(nx/20);
    patchLenBaseY= round(ny/20);
else
   patchLenBaseX= patchLenBase(1);
   patchLenBaseY= patchLenBase(2);
end

for ipatch = 1:nPatches
    %=====================================================================%
    % create arbitrary patch for each outputs
    %=====================================================================%   
    patchLenX = patchLenBaseX * (2+randsample(5,1)); % patchLen \in K*{4,5,6}
    patchLenY = patchLenBaseY * (2+randsample(5,1)); % patchLen \in K*{4,5,6}
    
    %=====================================================================%
    % ensure patches don't overlap
    % - assign patch at interval region (1:p/npatch) + offset
    % - so length of each patch is less than this interval width
    % (interval width = p/npatch)
    %=====================================================================%
    intervalWidthX = round(nx/nPatches);
    intervalWidthY = round(ny/nPatches);
    
    offsetX = intervalWidthX*(ipatch-1);
    offsetY = intervalWidthY*(ipatch-1);
    
    patchStartX = offsetX+randsample(intervalWidthX-patchLenX-1,1); %
    patchStartY = offsetY+randsample(intervalWidthY-patchLenY-1,1); %
    
    %         W(patchStart+1:patchStart+patchLen,iq) = tak_sample_signed_unif([5,10],1);
    tmp = zeros(nx,ny);
    tmp(patchStartX+1:patchStartX+patchLenX,patchStartY+1:patchStartY+patchLenY) = ...
        tak_sample_signed_unif([5,15],1);
%     keyboard
    W = W + tmp;
end
w=W(:);