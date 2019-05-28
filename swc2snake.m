function [ snakes ] = swc2snake( swc, I )
%SWC2SNAKE convert trace in SWC format to a population of snakes. One
%snakes per branch.

n = size(swc,1);
nb = swc(:,[1 7]);
nb = nb(all(nb>0, 2), :);

% convert to graph
dv = zeros(n,3);
dv(nb(:,1),:) = swc(nb(:,2),3:5) - swc(nb(:,1),3:5);
dv = dv ./ (max(sqrt(sum(dv.^2,2)), eps)*ones(1,3));
dv2sq = sum((dv(nb(:,1),:) - dv(nb(:,2),:)).^2, 2);
G = sparse(nb(:,1), nb(:,2), dv2sq+eps, n, n);

% find non-root leaf nodes
[~, dt, ft, pred] = dfs(G+G',1); 
leaves = find(dt+1 == ft & (1:n)'~=1);      

% sort path by branch curvature
branches = cell(length(leaves),1);
brlen = zeros(length(leaves),1);
for num = 1:length(leaves)
    leaf = leaves(num);
    branches{num} = path_from_pred(pred, leaf);
    pnt = [branches{num}(2:end); branches{num}(1:end-1)]';
    brlen(num) = full(sum( diag( G(pnt(:,1), pnt(:,2)) ) ));
end
[~, order] = sort(brlen);

% convert to snakes
visit = false(n,1);
snknum = zeros(n,1);
snakes = struct('vert', [], ...                 % points on snake
        'collide', 0,...                        % collided snk# & snk pnt#
        'num', num2cell(1:length(leaves)) ...   % snake number
        );
for num = 1:length(leaves)
    % find branch in swc file
    branch = branches{order(num)};
    st = find(visit(branch), 1, 'last');
    if isempty(branch(st)), 
        st = 1; 
    else
        snakes(num).collide = snknum(branch(st));
    end
    
    % remesh the branch so points are one pixel apart
    snakes(num).vert = swc(branch(st:end),3:5)+1;
    snakes(num) = AC_remesh(snakes(num));
    
    % if it is not root, remove the first point (overlap at intersection)
    if snakes(num).collide > 0
        snakes(num).vert = snakes(num).vert(2:end,:);
    end
    
    snknum(branch(st:end)) = num;
    
    visit(branch) = true;
end

% plot
if nargin == 2
    figure(1), imshow(I);
    hold on;
    plot(swc(1,3), swc(1,4), 'ro');
    for num = 1:length(leaves)
        plot(snakes(num).vert(:,1), snakes(num).vert(:,2), '-o', 'linewidth', 2);
        pause;
    end
    hold off;
end
end
