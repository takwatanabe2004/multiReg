function prox=tak_isoTV_prox_4d(x,b,tau)
% prox=tak_isoTV_prox_4d(x,b,tau)
% (02/13/2014)
%=========================================================================%
% tau <- prox scaling term
%=========================================================================%
%%
dim = 4;
p = length(x)/dim;

%=========================================================================%
% indices of the gradients for the 4 coordinate directions
%=========================================================================%
ix1=(1    :   p);
iy1=(1+  p: 2*p);
ix2=(1+2*p: 3*p);
iy2=(1+3*p: 4*p);

%=========================================================================%
% break the signal into 4 components 
% (each component represents a coordinate direction)
%=========================================================================%
dx1=x(ix1);
dy1=x(iy1);
dx2=x(ix2);
dy2=x(iy2);

% apply masking on the signal
bdx1=b(ix1).*dx1(:);
bdy1=b(iy1).*dy1(:);
bdx2=b(ix2).*dx2(:);
bdy2=b(iy2).*dy2(:);

%=========================================================================%
% signal norm with the appropriate masking applied: (p x 1)
%=========================================================================%
% vector of euclidean norm of the gradients, with the mask applied
xnorm_mask = sqrt(bdx1.^2 + bdy1.^2 + bdx2.^2 + bdy2.^2);

%=========================================================================%
% Compute shrinkage factor
%-------------------------------------------------------------------------%
% tmp = tau/||w|| for the shrinkage, with the convention (0/0) = 0
%=========================================================================%
tmp=tau./xnorm_mask;
tmp(isnan(tmp))=0;
shrink=max(1-tmp,0); % shrinkage factor (4*ptil x 1)

%=========================================================================%
% apply shrinkage
%=========================================================================%
bdx1=shrink.*bdx1;
bdy1=shrink.*bdy1;
bdx2=shrink.*bdx2;
bdy2=shrink.*bdy2;
shrinked_signal=cat(1,bdx1,bdy1,bdx2,bdy2);

%=========================================================================%
% final Prox output
%-------------------------------------------------------------------------%
% - the prox output has two parts: masked part and non-masked parts, where
%  the masked part of the prox simply returns the input signal
%-------------------------------------------------------------------------%
prox=x;

% the non-masked part of the prox output; here return the shrinked signal
prox(b)=shrinked_signal(b);