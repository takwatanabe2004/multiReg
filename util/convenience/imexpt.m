function imexpt
% 07/10/2013 -> on my eecs computer, i want the display to popout on my secondary monitor
host=tak_get_host;
if strcmpi(host,'takanori.eecs.umich.edu')
    % 07/10/2013 my eecs computer...display on my secondary monitor
    offy=90;
    tak_imset2(gcf,1921,2.7*offy + (900-offy)/2,1440,(860-offy)/2.15) % hp2760p
elseif strcmpi(host,'takanori-PC')
    set(gcf,'Units','pixels','Position', [1 570 1920 432])
else
    set(gcf,'Units','pixels','Position', [1 600 1920 432])
end