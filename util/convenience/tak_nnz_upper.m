function [nnz_upper, sp_level, n_entries] = tak_nnz_upper(A)
% nnz_upper = tak_nnz_upper(A)
%|--------------------------------------------------------------------------------------|%
%| Description:
%|      Count the number of nonzero entries in the upper triangular part
%|      of the matrix A (hence, not counting the diagonals)
%|
%|--------------------------------------------------------------------------------------|%
%| Input: A = square matrix
%|
%| Output: nnz_upper = # nonzeroes in the upper-triangular part of A
%|         sp_level = sparsity level, ie, the fraction of the entires
%|                    in the upper-triangular that were nonzero.
%|        n_entries = nchoosek( size(A,1),2)
%|                  => the # entries in the upper triangular part of A      
%|--------------------------------------------------------------------------------------|%
%| Created 8/08/2012
%|--------------------------------------------------------------------------------------|%
%%
nnz_upper = full(sum(sum(triu(A~=0,1))));

p = size(A,1);
n_entries = nchoosek(p,2);
sp_level = nnz_upper/n_entries;