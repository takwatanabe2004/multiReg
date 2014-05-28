% vector shrinkage
% 6/25/2013
function z = tak_softvec(x, tau)
    z = max(1 - tau/norm(x),0)*x;
end