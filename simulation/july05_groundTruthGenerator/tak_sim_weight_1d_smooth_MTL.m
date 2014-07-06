function [W,idx_supp] = tak_sim_weight_1d_smooth_MTL(p,q)
% [W,idx_supp] = tak_sim_weight_1d_smooth_MTL(p,q)
%=========================================================================%
% - Generate 1d grouth truth weight vector support.
%=========================================================================%
% (07/05/2014)
%%
 
%=========================================================================%
% set smooth sinusoidal pulse signal parameters
%=========================================================================%
% length of sinusoidal pulse
L=p/8;

% frequency
freq1 = 2/(1*L);
freq2 = 1/(1*L);
freq3 = 2/(1*L);

idx_support = cell(4,1); % 4 sinusoidal pulse
idx_support{1} = (1:L) + round(p*0.1);
idx_support{2} = (1:L) + round(p*0.35);
idx_support{3} = (1:L) + round(p*0.7);
% idx_support{4} = (1:k) + round(n*0.9);

%=========================================================================%
% assign weights
%=========================================================================%
W = zeros(p,q);
for k=1:q
    % magnitude of sinusoidal pulses
    mag = tak_sample_signed_unif([4,8],4);
%     keyboard
    W(idx_support{1},k) = mag(1)*sin(freq1*pi*(1:L))';
    W(idx_support{2},k) = mag(2)*sin(freq2*pi*(1:L))';
    W(idx_support{3},k) = mag(3)*sin(freq3*pi*(1:L))';
    % w(idx_support{4}) = mag(4)*sin(.02*pi*(1:k));
end
idx_supp = find(W(:,1));