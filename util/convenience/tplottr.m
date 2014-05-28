function tplottr(x,y,options)
% used to be tplot.m
% 6/17/2013 -> added xlim default
if(~exist('options','var')||isempty(options))
    options={'linewidth',3};
end
figure,imexpr
if nargin==1
    plot(x,options{:})
    xlim([1,length(x)])
else
    plot(x,y,options{:})
    xlim([min(x),max(x)])
end
title(inputname(1),'Interpreter','none')
grid on
drawnow