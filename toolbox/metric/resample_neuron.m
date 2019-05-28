%RESAMPLE_NEURON resampling neuron points so that there is a
%unit vector of distance 1 between each node
%Purpose: allows MATLAB to go through the data quickly and more efficiently
%while still calculating a reasonably accurate prediction of the data 
%points to be used for nodes
%Input: swc neuron file
%Output: interpolated neuron node points

function newdata = resample_neuron(data)
    newdata = data(1,:);
    % do BFS with any point as root node
    n = size(data,1);
    idx = data(:,7) > 0;
    nb = data(:,[1 7]);
    nb = nb(idx,:);
    G = sparse(nb(:,1), nb(:,2),1,n,n);
    G = G + G';
    
%     xx = [data(nb(:,1),3), data(nb(:,2),3), nan(size(nb,1),1)]';
%     yy = [data(nb(:,1),4), data(nb(:,2),4), nan(size(nb,1),1)]';
%     figure(1), subplot(1,2,1); plot(xx(:), yy(:), '-o', 'LineWidth', 4);
%     set(gca,'Ydir','reverse')

    [d, ~, pred] = bfs(G,1);
    [~,order] = sort(d);

    deg = histc(nb(:),1:n);
    visit = zeros(n,1); visit(1) = 1;
    % go through each segment (curve between branchpoints or branch and 
    % leaf) in order sorted by the distance from root (first point)
    for i = find(deg(order)~=2)'
        index = order(i);
        % find path from leaf or branch node to root
        path = path_from_pred(pred,index);
        t = find(visit(path), 1, 'last');
        if isempty(t)
            error('no previous point, impossible');
        end
        % trim path to segment
        path = path(t:end);
        if length(path) <= 1, continue; end;    %
        % interpolate it, so distance between nodes is one
        swcpath = data(path, :);
        v = sqrt(sum( (swcpath(1:end-1,3:5) - swcpath(2:end,3:5)).^2,2));
        swcpath(v==0,:) = [];
        v(v==0) = [];
        f = [0; cumsum(v)];
        N = ceil(sum(v))*2+3;         % number of sampling pnts in a segment
        xi = linspace(0,sum(v),N);
            
        points = zeros(N,3);
        for h = 1:3
            points(:,h) = interp1(f,swcpath(:,h+2),xi);
        end
        % put back to swc format
        np = size(points,1);
        [~, pos] = min(abs(bsxfun(@minus, f, xi)));
        newd = [(0:np-1)' swcpath(pos,2) points swcpath(pos,6) (-1:np-2)'];
        newd(1,:) = [];
        newd(:,[1 7]) = newd(:,[1 7]) + size(newdata,1);
        newd(1,7) = visit(path(1));
        newdata = [newdata; newd];

        visit(index) = size(newdata,1); 
    end
%     keyboard;
% 
%     idx = newdata(:,7) > 0;
%     nb = newdata(:,[1 7]);
%     nb = nb(idx,:);
%     xx = [newdata(nb(:,1),3), newdata(nb(:,2),3), nan(size(nb,1),1)]';
%     yy = [newdata(nb(:,1),4), newdata(nb(:,2),4), nan(size(nb,1),1)]';
%     subplot(1,2,2); plot(xx(:), yy(:), '-o', 'LineWidth', 4);
%     set(gca,'Ydir','reverse')
end