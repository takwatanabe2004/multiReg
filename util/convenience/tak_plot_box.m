function tak_plot_box(boxroi, linecolor, linewidth)

if(~exist('linecolor','var')||isempty(linecolor))
    linecolor = 'm'; 
end

if(~exist('linewidth','var')||isempty(linewidth))
    linewidth = 2; 
end

x1 = boxroi(1);
x2 = boxroi(2);
y1 = boxroi(3);
y2 = boxroi(4);

hold on
line( [x1,x1], [y1,y2], 'linewidth',linewidth, 'color',linecolor)
line( [x2,x2], [y1,y2], 'linewidth',linewidth, 'color',linecolor)
line( [x1,x2], [y1,y1], 'linewidth',linewidth, 'color',linecolor)
line( [x1,x2], [y2,y2], 'linewidth',linewidth, 'color',linecolor)
