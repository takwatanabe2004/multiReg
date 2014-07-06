function [W,idx_supp] = tak_sim_weight_1d_patch_MTL(p,q)
% [W,idx_supp] = tak_sim_weight_1d_patch_MTL(p,q)
%=========================================================================%
% - Generate 1d grouth truth weight vector support.
%=========================================================================%
% (07/05/2014)
%%
 
%=========================================================================%
% set patch paramters
%=========================================================================%
% patch lengths
patch1 = round(0.07*p);
patch2 = round(0.1*p);
patch3 = round(0.08*p);
patch4 = round(0.05*p);

% patch offeets
offsets1 = round(0.1*p);
offsets2 = round(0.3*p);
offsets3 = round(0.6*p);
offsets4 = round(0.8*p);

idx_support = cell(4,1); % 4 sinusoidal pulse
idx_support{1} = (1:patch1) + offsets1;
idx_support{2} = (1:patch2) + offsets2;
idx_support{3} = (1:patch3) + offsets3;
idx_support{4} = (1:patch4) + offsets4;

%=========================================================================%
% assign weights
%=========================================================================%
W = zeros(p,q);
for k=1:q
    % magnitude of tha patches
    mag = tak_sample_signed_unif([2,5],4);
    
    W(idx_support{1},k) = mag(1);
    W(idx_support{2},k) = mag(2);
    W(idx_support{3},k) = mag(3);
    W(idx_support{4},k) = mag(4);
end
idx_supp = find(W(:,1));