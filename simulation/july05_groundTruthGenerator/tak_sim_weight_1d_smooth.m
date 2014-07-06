function [w,idx_supp] = tak_sim_weight_1d_smooth(n)
% [w,idx_supp] = tak_sim_weight_1d_smooth(n)
%=========================================================================%
% - Generate 1d grouth truth weight vector support.
%=========================================================================%
% (07/05/2014)
%%
 
%=========================================================================%
% set smooth sinusoidal pulse signal parameters
%=========================================================================%
% length of sinusoidal pulse
k=n/8;

% frequency
freq1 = 2/(1*k);
freq2 = 1/(1*k);
freq3 = 2/(1*k);

idx_support = cell(4,1); % 4 sinusoidal pulse
idx_support{1} = (1:k) + round(n*0.1);
idx_support{2} = (1:k) + round(n*0.35);
idx_support{3} = (1:k) + round(n*0.7);
% idx_support{4} = (1:k) + round(n*0.9);

% magnitude of sinusoidal pulses
mag = tak_sample_signed_unif([4,8],4);

%=========================================================================%
% assign weights
%=========================================================================%
w = zeros(n,1);
w(idx_support{1}) = mag(1)*sin(freq1*pi*(1:k));
w(idx_support{2}) = mag(2)*sin(freq2*pi*(1:k));
w(idx_support{3}) = mag(3)*sin(freq3*pi*(1:k));
% w(idx_support{4}) = mag(4)*sin(.02*pi*(1:k));

idx_supp = find(w);