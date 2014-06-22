function [y dy] = tak_bspline4(x)
% I evaluated the 6 signal components using the equation from pdf-pg181 on 556 book.
% I'm sure I can evaluate a more compact expression using the absolute value function, but for
% now, I implement the brute force summation edition as in here.
% 12/13/2011 - Takanori, Research Noe 12/13/2011 pg.11

s1 =     (x + 5/2).^4    .*  (x >= -5/2);     
s2 =  -5*(x + 3/2).^4    .*  (x >= -3/2);      
s3 =  10*(x + 1/2).^4    .*  (x >= -1/2);      
s4 = -10*(x - 1/2).^4    .*  (x >=  1/2);      
s5 =   5*(x - 3/2).^4    .*  (x >=  3/2);      
s6 =    -(x - 5/2).^4    .*  (x >=  5/2);      

y = (s1 + s2 + s3 + s4 + s5 + s6)/24;

if nargout == 2
    dy = tak_bspline3(x + 1/2) - tak_bspline3(x - 1/2);
end
