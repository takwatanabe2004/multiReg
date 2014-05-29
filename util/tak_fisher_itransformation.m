function r = tak_fisher_itransformation(z)
%| inverse fisher transformation
r = (exp(2*z)-1)./(exp(2*z)+1);