function [y dy] = tak_bspline1(x)
y = [1-abs(x)].*(abs(x)<1);

if nargout == 2
    dy = tak_bspline0(x + 1/2) - tak_bspline0(x - 1/2);
end
