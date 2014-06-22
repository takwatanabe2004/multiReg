function [y dy] = tak_bspline3(x)
y = [2/3-x.^2+abs(x).^3/2]   .*   [abs(x) < 1] + ...
    [(2-abs(x)).^3/6]   .*   [1 <= abs(x) & abs(x) < 2];

if nargout == 2
    dy = tak_bspline2(x + 1/2) - tak_bspline2(x - 1/2);
end

% function out = tbsp(x)
% if abs(x) <1
%     out = 2/3-x^2+abs(x)^3/2;
% elseif abs(x)>=1 && abs(x)<2
%     out = (2-abs(x))^3/6;
% else
%     out = 0;
% end
% 
% function y = tak_bspline2(x)
% y = (3/4 - x.^2) .* (abs(x) < 1/2) + ...
% 	(abs(x) - 3/2).^2 / 2 .* (1/2 <= abs(x) & abs(x) < 1.5);
% 
% if abs(x) < 1/2
%     out = 3/4-x^2;
% elseif abs(x)>=1/2 && abs(x)<3/2
%     out = x^2/2-3/2*abs(x)+9/8;
% else
%     out = 0;
% end
