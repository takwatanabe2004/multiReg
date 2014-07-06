function [w,idx_supp] = tak_sim_weight_1d_patch(n)
% [w,idx_supp] = tak_sim_weight_1d_patch(n)
%=========================================================================%
% - Generate 1d grouth truth weight vector support.
%=========================================================================%
% (07/05/2014)
%%
 
%=========================================================================%
% set patch paramters
%=========================================================================%
% patch lengths
patch1 = round(0.07*n);
patch2 = round(0.1*n);
patch3 = round(0.08*n);
patch4 = round(0.05*n);

% patch offeets
offsets1 = round(0.1*n);
offsets2 = round(0.3*n);
offsets3 = round(0.6*n);
offsets4 = round(0.8*n);

idx_support = cell(4,1); % 4 sinusoidal pulse
idx_support{1} = (1:patch1) + offsets1;
idx_support{2} = (1:patch2) + offsets2;
idx_support{3} = (1:patch3) + offsets3;
idx_support{4} = (1:patch4) + offsets4;

% magnitude of tha patches
mag = tak_sample_signed_unif([4,8],4);

%=========================================================================%
% assign weights
%=========================================================================%
w = zeros(n,1);
w(idx_support{1}) = mag(1);
w(idx_support{2}) = mag(2);
w(idx_support{3}) = mag(3);
w(idx_support{4}) = mag(4);

idx_supp = find(w);