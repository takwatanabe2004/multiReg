function tstemm2(x1,x2,options)
% for superimposed stem-plot of two signals of same length
%----------------------------------------------------------------------------------
% The function can be called in the following ways:
% tplot(x)
% tplot(x,y)
% tplot(x,options)
% tplot(x,y,options)
%----------------------------------------------------------------------------------
% (06/09/2014)
%%
figure
if nargin<3
    options={'linewidth',2}; % <-default
end
stem(x1,options{:},'color','b');
hold on
stem(x2,options{:},'color','r');
xlim([1,length(x1)])

title(inputname(1),'Interpreter','none')
grid on
drawnow