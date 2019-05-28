function [e, l] = fusion_moves (U, Ps, edge_ends, max_num_iterations, verbose)
% Perform alpha-expansion moves using QPBO, also known as fusion moves.
% Parameters
% ----------
% U : array of shape = (num_nodes, num_labels)
%     The matrix of unary energies, where `U(i, j)` is the energy of node
%     `i` taking label `j`.
% Ps : array of shape = (num_labels, num_labels, num_edges)
%     The list of matrices of pairwise energies, where `P(i, j, k)` is the
%     energy of adjacent nodes taking labels `i` and `j`, according to the
%     edge `k`.
% edge_ends : array of shape = (num_edges, 2)
%     where row `edge_ends(k,:)` is the indices of two endpoints of the
%     edge k
% max_num_iterations : int
%     The maximum number of outer iterations.
% verbose : binary
%     true for showing message and false to display nothing
% Returns
% -------
% e : float
%     The final energy.
% l : array of shape = (num_nodes,1)
%     The array of refined labels corresponding to `e`.

if nargin < 3
    error('Need more inputs');
end
if nargin < 4
    max_num_iterations = 10;
end
if nargin < 5
    verbose = false;
end

%% input validation
num_nodes = size(U,1);
num_labels = size(U,2);
num_edges = size(edge_ends,1);

if ~isa(edge_ends, 'double')
    edge_ends = double(edge_ends);
end

if size(edge_ends,2) ~= 2,
    error('edge_ends must have size [num_edges, 2]');
end
if size(Ps,1) ~= size(Ps,2) || size(Ps,1) ~= num_labels || size(Ps,3) ~= num_edges
    error(['Ps must have size [num_labels, num_labels, num_edges] = [ ' ...
            num2str([num_labels num_labels num_edges]) ' ] ']);
end

if max_num_iterations <= 0
    error('max_num_iterations must be greater than zero');
end

%% perform fusion moves
% init label
l = ceil(num_labels/2)*ones(num_nodes,1);
% init energy
last_e = get_energy(U, Ps, edge_ends, l);
% init unary and pairwise terms
terminalWeights = zeros(num_nodes,2);
edgeWeights = zeros(num_edges,6);
edgeWeights(:,1:2) = edge_ends;
% start fusion moves
for iter = 1:max_num_iterations
    % For each `alpha`, solve a binary sub-problem. The binary sub-problem is:
    % Do non-`alpha` nodes keep their current label or change to `alpha`?
    for alpha = 1:num_labels        
        % `E0` is the energy of retaining the current label;
        % `E1` is the energy of node `j` changing to label `alpha`.
        E0idx = (1:num_nodes)' + (l-1)*num_nodes;
        
        % `E00` = P(l[k],  l[j]),  `E01` = P(l[k],  alpha), 
        % `E10` = P(alpha, l[j]),  `E11` = P(alpha, alpha));
        E00idx = l(edge_ends(:,1)) + (l(edge_ends(:,2))-1)*num_labels + (0:num_edges-1)'*num_labels^2;%sub2ind(size(Ps), l(edge_ends(:,1)), l(edge_ends(:,2)), (1:num_edges)');
        E01idx = l(edge_ends(:,1)) + (alpha-1)*num_labels + (0:num_edges-1)'*num_labels^2;%sub2ind(size(Ps), l(edge_ends(:,1)), alpha*ones(num_edges,1), (1:num_edges)');
        E10idx = alpha + (l(edge_ends(:,2))-1)*num_labels + (0:num_edges-1)'*num_labels^2;%sub2ind(size(Ps), alpha*ones(num_edges,1), l(edge_ends(:,2)), (1:num_edges)');
        
        % apply QPBO
        terminalWeights(:,1) = U(E0idx);        % E0
        terminalWeights(:,2) = U(:,alpha);      % E1
        
        edgeWeights(:,3) = Ps(E00idx);          % E00
        edgeWeights(:,4) = Ps(E01idx);          % E01
        edgeWeights(:,5) = Ps(E10idx);          % E10
        edgeWeights(:,6) = Ps(alpha,alpha,:);   % E11

        [~, labels] = qpboMex(terminalWeights, edgeWeights);
        % switch label
        l(labels == 1) = alpha;
    end

    e = get_energy(U, Ps, edge_ends, l);
    % check if converge
    if e >= last_e, break; end
    last_e = e;
end

if verbose
    if e < last_e
        disp(['Did not converge at iter ' num2str(iter)]);
    else
        disp(['Converge at iter ' num2str(iter)]);
    end
end
    %% auxiliarly functions
    function e = get_energy(U, Ps, edge_ends, l)
        uidx = (1:size(U,1))' + (l-1)*size(U,1);
        eidx = l(edge_ends(:,1)) + (l(edge_ends(:,2))-1)*size(Ps,1) + (0:size(edge_ends,1)-1)'*size(Ps,1)*size(Ps,2);
        e = sum(U(uidx)) + sum(Ps(eidx));
    end
end