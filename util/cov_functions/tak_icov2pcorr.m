function [pcorr,pcorr2,pcorr3,pcorr4] = tak_icov2pcorr(icov)
% pcor = tak_icov2pcorr(icov)
%|------------------------------------------------------------------------------|%
%| Input: 
%|   icov = inverse-covariance matrix
%|------------------------------------------------------------------------------|%
%| Output: 
%|   pcorr = partial-correlation matrix (think of it as conditional correlation)
%|------------------------------------------------------------------------------|%
%|  Convert the inverse-covariance into partial correlation.
%|  Current implementation is rather crude, but in GGM we focus on sparse
%|  icov model, so actual matrix-multiplication computation is not too intensive. 
%|------------------------------------------------------------------------------|%
%| 9/15/2012
%| 04/14/2013 -> error???.  see pg.6, equation (7) of the following:  http://arxiv.org/pdf/1107.1270v3.pdf
%|              (though they set the diagonals of the pcorr matrix to zero, whereas we set it to 1 here)
%| 05/13/2013 -> freewill does not have 'cov2corr' from finance toolbox...
%|               saved an mcentral script for it
%%
D = diag((diag(icov)).^(-1/2));
pcorr = D*icov*D;

% use Cov2Corr from mcentral
[trash,pcorr2] = Cov2Corr(icov);

p = size(icov,1);
% flip the signs in the off-diagonals
mask_offdiag = ~logical(eye(p));

pcorr3 = pcorr2;
pcorr3(mask_offdiag)=-pcorr3(mask_offdiag);

pcorr4 = zeros(p);
for i = 1:p
    for j = 1:p
        if i==j
            pcorr4(i,j) = 1;
        else
            pcorr4(i,j) = -icov(i,j)/sqrt(icov(i,i)*icov(j,j));
        end
    end
end