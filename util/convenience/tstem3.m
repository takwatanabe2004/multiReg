function tstem3(x1,x2,x3,options)
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
if nargin<4
    options={'linewidth',2}; % <-default
end
stem(x1,options{:},'color','b');
hold on
stem(x2,options{:},'color','r');
stem(x3,options{:},'color',[0 0.8 0]);
xlim([1,max([length(x1),length(x2),length(x3)])])

grid on
drawnow