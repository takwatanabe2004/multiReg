function [nnz_lower, sp_level, n_entries] = tak_nnz_lower(A)
% nnz_lower = tak_nnz_lower(A)
%|--------------------------------------------------------------------------------------|%
%| Description:
%|      Count the number of nonzero entries in the lower triangular part
%|      of the matrix A (hence, not counting the diagonals)
%|
%|--------------------------------------------------------------------------------------|%
%| Input: A = square matrix
%|
%| Output: nnz_upper = # nonzeroes in the lower-triangular part of A
%|         sp_level = sparsity level, ie, the fraction of the entires
%|                    in the lower-triangular that were nonzero.
%|        n_entries = nchoosek( size(A,1),2)
%|                  => the # entries in the lower triangular part of A      
%|--------------------------------------------------------------------------------------|%
%| Created (06/11/2014)
%|--------------------------------------------------------------------------------------|%
%%
nnz_lower = full(sum(sum(tril(A~=0,-1))));

p = size(A,1);
n_entries = nchoosek(p,2);
sp_level = nnz_lower/n_entries;