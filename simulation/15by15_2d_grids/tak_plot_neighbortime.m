function tak_plot_neighbortime(SIGNAL,x,y,YLIM)
% (06/27/2014)
% Brought this from my old "structSVM" repos.
% One major change: SIGNAL = (X, Y, T)...before it was (T,X,Y)
%%
% (12/22/2013)
% - created so that i can plot the nearest neighbor timeseries of the
%   spatiotemporal signal i create
%%
% SIGNAL = (T, nx, ny) spatiotemporal signal sampled from matrix normal 
%           distribution, specifically, created from tak_prec_connectome2d.m
% (x,y) = the location you want to plot the nearest neighbor timeseries
% YLIM = (optional) ylimit axis
%%
figure,imexpb
subplot(211)
plot(squeeze(SIGNAL(x,y,:)),'linewidth',4),hold on,
xlim([1,size(SIGNAL,3)]);
% AXIS=axis;
plot(squeeze(SIGNAL(x-1,y,:)),'r','linewidth',4)
plot(squeeze(SIGNAL(x+1,y,:)),'g','linewidth',4)
legend('CENTER','x-1','x+1')
title('(x-1,y)...(x+1,y)...')
grid on

subplot(212)
plot(squeeze(SIGNAL(x,y,:)),'linewidth',4),hold on
plot(squeeze(SIGNAL(x,y+1,:)),'k','linewidth',4)
plot(squeeze(SIGNAL(x,y-1,:)),'m','linewidth',4)
legend('CENTER','y-1','y+1')
title('(x,y-1)...(x,y+1)...')
grid on
xlim([1,size(SIGNAL,3)]);

if nargin == 4 && ~isempty(YLIM)
    subplot(211),ylim(YLIM)
    subplot(212),ylim(YLIM)
end
% keyboard