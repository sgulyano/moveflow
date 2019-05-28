function data = graph2swc( G, N )
%GRAPH2SWC convert graph as sparse adjency matrix and nodes (xyz + radius)
%into SWC format. Assume G is one connected tree 
[~,dt,pred] = bfs(G,1);
[b,ix] = sort(dt);
order = ix(b>0)';

nc = size(N,1);
data = [(1:nc)', 2*ones(nc,1), N, pred'];
data = data(order,:);
% get fixed index
fixind = zeros(size(order));
fixind(order) = 1:length(order);
% change to fixed index
data(:,1) = fixind(data(:,1));
data(1,7) = -1;
data(2:end,7) = fixind(data(2:end,7));
end

