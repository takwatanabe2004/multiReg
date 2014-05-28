function imexptr
% 07/10/2013 -> on my eecs computer, i want the display to popout on my secondary monitor
% set(gcf,'Units','pixels','Position', [1+1920/2 1+1080/2 1920/2 1080/2])
% set(gcf,'Units','pixels','Position', [1+1920/2 600 950 432])

% 07/10/2013 -> on my eecs computer, i want the display to popout on my secondary monitor
% set(gcf,'Units','pixels','Position', [1 1+1080/2 1920/2 1080/2])
% set(gcf,'Units','pixels','Position', [1 600 950 432])

host=tak_get_host;
if strcmpi(host,'takanori.eecs.umich.edu')
    % 07/10/2013 my eecs computer...display on my secondary monitor
    offy=90;
    tak_imset2(gcf,1921+1440/2,2.7*offy + (900-offy)/2,1440/2,(860-offy)/2.15) % hp2760p
else
    set(gcf,'Units','pixels','Position', [1+1920/2 568 950 440]) % for my windows7 lab computer
end

