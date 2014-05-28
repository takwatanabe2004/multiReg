function tplot(x,varargin)
% the old tplot.m is now called tplott.m
%----------------------------------------------------------------------------------
% The function can be called in the following ways:
% tplot(x)
% tplot(x,y)
% tplot(x,options)
% tplot(x,y,options)
%----------------------------------------------------------------------------------
% 6/15/2013
% 6/17/2013 -> added xlim default
% 6/24/2013 -> began using varargin
%%
options={'linewidth',2}; % <-default
switch length(varargin)
    case 0
        plot(x,options{:});
        xlim([1,length(x)])
    case 1
        if isa(varargin{1},'cell')
            options=varargin{1};
            plot(x,options{:})
            xlim([1,length(x)])
        else
            y=varargin{1};
            plot(x,y,options{:})
            xlim([min(x),max(x)])
        end
    case 2
        y=varargin{1};
        options=varargin{2};
        plot(x,y,options{:})
        xlim([min(x),max(x)])
    otherwise
        error('invalid input set')
end

title(inputname(1),'Interpreter','none')
grid on
drawnow