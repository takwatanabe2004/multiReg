function [y, dy] = tak_bspline5(x)
% I evaluated the 7 signal components using the equation from pdf-pg181 on 556 book.
% I'm sure I can evaluate a more compact expression using the absolute value function, but for
% now, I implement the brute force summation edition as in here.
% 12/13/2011 - Takanori, Research Noe 12/13/2011 pg.11

s1 =     (x + 3).^5    .*  (x >= -3);     
s2 =  -6*(x + 2).^5    .*  (x >= -2);      
s3 =  15*(x + 1).^5    .*  (x >= -1);      
s4 = -20*(  x  ).^5    .*  (x >=  0);      
s5 =  15*(x - 1).^5    .*  (x >=  1);      
s6 =  -6*(x - 2).^5    .*  (x >=  2);      
s7 =     (x - 3).^5    .*  (x >=  3);      

y = (s1 + s2 + s3 + s4 + s5 + s6 + s7)/120;

if nargout == 2
    dy = tak_bspline4(x + 1/2) - tak_bspline4(x - 1/2);
end
