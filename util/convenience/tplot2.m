function tplot2(x1,x2,options)
% for superimposed plot of two signals of same length
%----------------------------------------------------------------------------------
% The function can be called in the following ways:
% tplot(x)
% tplot(x,y)
% tplot(x,options)
% tplot(x,y,options)
%----------------------------------------------------------------------------------
% (05/28/2014)
% (06/18/2014) - modification on xlim
%%
if nargin<3
    options={'linewidth',2}; % <-default
end
plot(x1,options{:},'color','b');
hold on
plot(x2,options{:},'color','r');
% xlim([1,length(x1)])
xlim([1,max(length(x1),length(x2))]) % (06/18/2014) modification

title(inputname(1),'Interpreter','none')
grid on
drawnow