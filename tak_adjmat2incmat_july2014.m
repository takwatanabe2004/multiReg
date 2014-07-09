function C=tak_adjmat2incmat_july2014(A,option)
%=========================================================================%
% - modified so that C contains the signed correlation info from the
%   weighted-adjacency matrix A
%-------------------------------------------------------------------------%
% option = {0,1,2}
%   - 0: binary   signed weight (default)
%   - 1: absolute signed weight
%   - 2: squared  signed weight
%=========================================================================%
% C=tak_adjmat2incmat(A)
% - convert from |V|x|V| adjacency matrix A to |E|x|V| incidence matrix.
% - V = nodes, E = edges.
% (code from http://www.mathworks.com/matlabcentral/fileexchange/24661-graph-adjacency-matrix-to-incidence-matrix/content/adj2inc.m)
%=========================================================================%
% (07/08/2014)
%%
if nargin==1
    option=0;
end
n_nodes=size(A,1);
[vNodes1,vNodes2] = find(triu(A~=0));     
n_edges=length(vNodes1);
% keyboard

% Recall every row ni the incidence matrix contains a -1 and a +1 entry.
% first  half takes care of the {-1} entry in the edge graph rows
% second half takes care of the {+1} entry in the edge graph rows
term1=[1:n_edges, 1:n_edges]';
term2=[vNodes1; vNodes2];

% term3=[-ones(n_edges,1);ones(n_edges,1)];
term3 = zeros(2*n_edges,1);
% term3(n_edges+1:end)=1;
for ii=1:n_edges
    rval = A(vNodes1(ii),vNodes2(ii));
    
    switch option
        case 0
        %=================================================================%
        % unweighted case
        %=================================================================%
        term3(ii) = -sign(rval);
        term3(ii+n_edges) = 1;
    
        case 1
        %=================================================================%
        % absolute vaule weighted case
        %=================================================================%    
        term3(ii) = -sign(rval)*abs(rval);
        term3(ii+n_edges) = abs(rval);

        case 2
        %=================================================================%
        % squared vaule weighted case
        %=================================================================%
        rval = A(vNodes1(ii),vNodes2(ii));
        term3(ii) = -sign(rval)*(rval^2);
        term3(ii+n_edges) = (rval)^2;
    end
end

C = sparse(term1, term2, term3, n_edges, n_nodes);
