function [w,W,idx_supp] = tak_sim_weight_2d_patch(nx,ny)
% [w,W,idx_supp] = tak_sim_weight_2d_patch(n)
%=========================================================================%
% - Generate 2d grouth truth weight vector.
%=========================================================================%
% (07/05/2014)
%% 

%=========================================================================%
% set patch signal parameters
%=========================================================================%
% size of each patches
nsize1 = round([0.1*nx 0.2*ny]);
nsize2 = round([0.25*nx 0.25*ny]);
nsize3 = round([0.2*nx 0.1*ny]);
% nsize{4} = round([10 10]);

% offsets locations of patches
offset1 = round([nx*.2, ny*.5]);
offset2 = round([nx*.5, ny*.1]);
offset3 = round([nx*.6, ny*.7]);

% indices of the patches
xidx1 = offset1(1) + (1:nsize1(1));   yidx1 = offset1(2) + (1:nsize1(2));
xidx2 = offset2(1) + (1:nsize2(1));   yidx2 = offset2(2) + (1:nsize2(2));
xidx3 = offset3(1) + (1:nsize3(1));   yidx3 = offset3(2) + (1:nsize3(2));

% magnitude of sinusoidal pulses
mag = tak_sample_signed_unif([4,8],4);
% keyboard
%=========================================================================%
% assign weights
%=========================================================================%
W = zeros(nx,ny);
W(xidx1,yidx1) = mag(1);
W(xidx2,yidx2) = mag(2);
W(xidx3,yidx3) = mag(3);
% w(idx_support{4}) = mag(4)*sin(.02*pi*(1:k));

w=W(:);
idx_supp = find(w);