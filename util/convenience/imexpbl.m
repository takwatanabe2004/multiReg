function imexpbl
% 07/10/2013 -> on my eecs computer, i want the display to popout on my secondary monitor
% set(gcf,'Units','pixels','Position', [1 1 1920/2 1080/2])
% set(gcf,'Units','pixels','Position', [1 1 950 432])

host = tak_get_host;

if strcmpi(host,'takanori-HP')
    tak_imset2(gcf,641,50,1280/2,670) % hp2760p
elseif strcmpi(host,'takanori.eecs.umich.edu')
    % 07/10/2013 my eecs computer...display on my secondary monitor
    offy=90;
    tak_imset2(gcf,1921,offy,1440/2,(860-offy)/2.15) % hp2760p
else
    set(gcf,'Units','pixels','Position', [1 41 950 440]) % for my windows7 lab computer
end
