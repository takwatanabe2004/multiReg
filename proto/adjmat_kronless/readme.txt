(06/28/2014)
- try to make "kron"-less version of adjmat generator
(the kron operator is the computational bottleneck for now).
%=========================================================================%
- bad news....just by completing the 1-d case, i can already see this 
  "non-kron" approach ain't helping....the "sparse" function is the bottleneck,
   when the index set is long enough....this is probably why kron was taking 
   so long too....fml