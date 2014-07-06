function samp = tak_sample_signed_unif(range,nsamp)
% samp = tak_sample_signed_unif(range,nsamp)
%=========================================================================%
% - sample from: unif([-range(2),-range(1)} || [range(1),range(2)])
%=========================================================================%
% (07/05/2014)
%%
% ensure range(1)<range(2)
range = sort(range,'ascend');

bias=range(1);        
mag=range(2)-range(1);
samp = (mag*rand(nsamp,1)+bias) .* sign(randn(nsamp,1));