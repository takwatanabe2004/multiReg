function [P, A] = tak_prec_ar1d(N, r)
% [P, A] = tak_make_prec1d(N, r)
%|-----------------------------------------------------------------------|%
%| Description:
%|      Create the precision matrix corresponding to the causal 1d AR(1) model.
%|-----------------------------------------------------------------------|%
%| Input: 
%|      N = vector length
%|      r = smoothness parameter
%|-----------------------------------------------------------------------|%
%| Output: 
%|      P = precision matrix in 'sparse' format.  
%|          Do full(P\speye(N)) to recover covariance matrix.
%|      A = A'*A=P
%|-----------------------------------------------------------------------|%
%| http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
%| To sample from matrix normal distribution, simply do this
% z=randn(N,1);
% x=A\z;
%
%| self note: (12/22/2013)
%|      COV=(A^-1 * A^-T)=(A^T A)^-1
%|       => A^-1 *z ~ COV distribution, so do A\z to sample from AR1D distribution.
%|                    as we need B*B^T = COV, where B:=A^-1
%|-----------------------------------------------------------------------|%
%| Created 2/13/2012
%|-----------------------------------------------------------------------|%
%| ar-model has this scaling component
scale = sqrt(1-r^2);

%| index for the main diagonal entries
ind_col_main = 1:N;
ind_row_main = 1:N;

%| entries for the main diagonal
entries_main = ones(N,1);
%| this modification necessary for the zero edge condition
entries_main(1) = scale;

%| index for the lower subdiagonal entries
ind_col_subl = 1:N-1;
ind_row_subl = 2:N;

%| entries for the lower subdiagonal
entries_subl = -r*ones(N-1,1);

%| concatenate the indices and entries
ind_col = [ind_col_main(:); ind_col_subl(:)];
ind_row = [ind_row_main(:); ind_row_subl(:)];
entries = [entries_main(:); entries_subl(:)];


A = sparse(ind_row, ind_col, entries/scale);
P = A'*A;
