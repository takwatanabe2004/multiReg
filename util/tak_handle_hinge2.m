function hinge2 = tak_handle_hinge2
% function handle for squared hinge loss function
% - fields contains function value, gradient value, and proximal operator
%%
hinge2.func=@(t) max(0,1-t).^2;      
hinge2.grad=@(t) (t<=1).*(2*t-2);  % derivative
hinge2.prox=@(t,tau) (t<=1).*(t+2*tau)/(1+2*tau) + (t>1).*t;