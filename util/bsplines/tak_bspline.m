function [y dy]= tak_bspline(x, deg)
% 12/13/2011: allows user to specify the bspline order.  default cubic.
if nargin == 1
    deg = 3;
end

if deg == 1
    [y,dy] = tak_bspline1(x);
elseif deg == 2
    [y,dy] = tak_bspline2(x);
elseif deg == 3
    [y,dy] = tak_bspline3(x);
elseif deg == 4
    [y,dy] = tak_bspline4(x);
elseif deg ==5
    [y,dy] = tak_bspline5(x);
else
    disp('error: only degree 1 to 5 supported')
end
