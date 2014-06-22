function [y dy]= tak_bspline2(x)
y = (3/4 - x.^2) .* (abs(x) < 1/2) + ...
	(abs(x) - 3/2).^2 / 2 .* (1/2 <= abs(x) & abs(x) < 1.5);

if nargout == 2
    dy = tak_bspline1(x + 1/2) - tak_bspline1(x - 1/2);
end
