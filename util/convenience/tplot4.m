function tplot4(x1,x2,x3,x4,options)
% for superimposed plot of two signals of same length
%----------------------------------------------------------------------------------
% The function can be called in the following ways:
% tplot(x)
% tplot(x,y)
% tplot(x,options)
% tplot(x,y,options)
%----------------------------------------------------------------------------------
% (06/18/2014)
%%
if nargin<5
    options={'linewidth',2}; % <-default
end
plot(x1,options{:},'color','b');
hold on
plot(x2,options{:},'color','r');
plot(x3,options{:},'color',[0 0.8 0]);
plot(x4,options{:},'color','c');
xlim([1,max([length(x1),length(x2),length(x3),length(x4)])])

grid on
drawnow