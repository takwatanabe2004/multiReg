function imexp
% 5/8/2013 -> hostname option
% 07/10/2013 -> on my eecs computer, i want the display to popout on my secondary monitor
%%
% set(gcf,'Units','pixels','Position', get(0,'ScreenSize'))
% 10/20/2011: Modified so there are some margin on the bottom of the screen to see the command window.
%%
host = tak_get_host;

if strcmpi(host,'takanori-HP')
    tak_imset2(gcf,1,80,1280,640) % hp2760p
elseif strcmpi(host,'takanori.eecs.umich.edu')
    % 07/10/2013 my eecs computer...display on my secondary monitor
    offy=90;
    tak_imset2(gcf,1921,offy,1440,900-offy) % hp2760p
else
    set(gcf,'Units','pixels','Position', [1 121 1920 850]) % for my windows7 lab computer
end