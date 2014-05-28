function hinge1 = tak_handle_hinge1
% function handle for hinge loss function
% - fields contains function value, gradient value, and proximal operator
%%
hinge1.func=@(t) max(0,1-t);
hinge1.grad=@(t) (t<1).*(-1); % undefined at t=1
hinge1.prox=@(t,tau) (t<(1-tau)).*(t+tau) + ((1-tau)<=t).*(t<=1).*1 + (1<t).*t;