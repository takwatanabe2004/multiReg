function tstemm(x,varargin)
% stem-plot version of tplot
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
options={'linewidth',2}; % <-default
switch length(varargin)
    case 0
        stem(x,options{:});
        xlim([1,length(x)])
    case 1
        if isa(varargin{1},'cell')
            options=varargin{1};
            stem(x,options{:})
            xlim([1,length(x)])
        else
            y=varargin{1};
            stem(x,y,options{:})
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