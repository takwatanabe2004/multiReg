function tak_plot_neighbortime8nn(SIGNAL,x,y,YLIM)
% (06/27/2014)
% Same as tak_plot_neighbortime, but plot the timeseries on a single plot,
% with the timeseries from the 8-nearest-negibhro plotted
%%
figure,imexpb
plot(squeeze(SIGNAL(x,y,:)),'linewidth',6),hold on,
xlim([1,size(SIGNAL,3)]);
% AXIS=axis;
plot(squeeze(SIGNAL(x-1,y,:)),'r','linewidth',4)
plot(squeeze(SIGNAL(x+1,y,:)),'g','linewidth',4)
plot(squeeze(SIGNAL(x,y-1,:)),'k','linewidth',4)
plot(squeeze(SIGNAL(x,y+1,:)),'m','linewidth',4)
plot(squeeze(SIGNAL(x-1,y-1,:)),'y','linewidth',4)
plot(squeeze(SIGNAL(x+1,y-1,:)),'c','linewidth',4)
plot(squeeze(SIGNAL(x-1,y+1,:)),'color',[0.4,0,0.3],'linewidth',4)
plot(squeeze(SIGNAL(x+1,y+1,:)),'color',.7*[1,1,1],'linewidth',4)
legend('CENTER','x-1','x+1','y-1','y+1',...
                'x-1,y-1','x+1,y-1','x-1,y+1','x+1,y+1')
grid on
grid on
xlim([1,size(SIGNAL,3)]);

if nargin == 4 && ~isempty(YLIM)
    ylim(YLIM)
end
% keyboard