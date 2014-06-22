%=========================================================================%
% quick sanity check that C'C+I does not have a sparse inverse...
% (kinda like how the precisionMatrix of COV_AR is sparse, but its inverse is not)
%=========================================================================%
clear
p = 20;
C = tak_diffmat_1d(p,0);
L=full(C'*C+speye(p));
Linv=inv(L);
% imcovv(L)
% imcovv(Linv)

%-------------------------------------------------------------------------%
% like the AR covariance, we have a constant decay rate
%-------------------------------------------------------------------------%
Linv(2,1)/Linv(1,1)
Linv(3,1)/Linv(1,1)
Linv(4,1)/Linv(1,1)
'-------'
Linv(3,1)/Linv(2,1)
Linv(4,1)/Linv(2,1)
Linv(5,1)/Linv(2,1)