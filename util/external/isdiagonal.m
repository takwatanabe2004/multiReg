function result = isdiagonal (A)
%% 06/23/2013
% REF: http://www.mathworks.com/matlabcentral/newsreader/view_thread/299383
% help block (yada, yada, yada)
[I,J] = find(A);
if ~isempty(I)
  result = all(I == J);
else
  % make the simple choice that an all zero matrix
  % or an empty array is by definition diagonal, since
  % it has no non-zero off diagonals.
  result = true;
end