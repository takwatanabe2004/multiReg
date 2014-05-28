function [P, A] = tak_prec_ar4d( imdim, smoothness)
% [P, A] = tak_prec_ar4d(imdim,smoothness)
%|------------------------------------------------------------------------------|%
%| Description:
%|      Create the precision matrix corresponding to the causal 4d AR(1) model.
%|------------------------------------------------------------------------------|%
%| Input: 
%|      imdim = [n1,n2,n3,n4]
%|      smoothness = [r1,r2,r3,r4]
%|------------------------------------------------------------------------------|%
%| Created 6/10/2013
%|------------------------------------------------------------------------------|%
n1 = imdim(1);
n2 = imdim(2);
n3 = imdim(3);
n4 = imdim(4);

r1 = smoothness(1);
r2 = smoothness(2);
r3 = smoothness(3);
r4 = smoothness(4);
if nargout == 1
    P1 = tak_prec_ar1d(n1,r1);
    P2 = tak_prec_ar1d(n2,r2);
    P3 = tak_prec_ar1d(n3,r3);
    P4 = tak_prec_ar1d(n4,r4);
    P = kron(P4,kron(P3, kron(P2,P1)));
else
    [P1, A1] = tak_prec_ar1d(n1,r1);
    [P2, A2] = tak_prec_ar1d(n2,r2);
    [P3, A3] = tak_prec_ar1d(n3,r3);
    [P4, A4] = tak_prec_ar1d(n4,r4);
    P = kron(P4,kron(P3, kron(P2,P1)));
    A = kron(A4,kron(A3, kron(A2,A1)));
end
