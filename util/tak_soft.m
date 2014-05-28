function y = tak_soft(t, tau)
% soft-threshold operator
%%
% if sum(abs(tau(:)))==0
%    y = x;
% else
%    y = max(abs(x) - tau, 0);
%    y = y./(y+tau) .* x;
% end
%%
% y = t.*max(0,1-tau./abs(t));
%%
y = sign(t).*max(0,abs(t)-tau);