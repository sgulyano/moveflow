function [ swc ] = snake2swc( snks )
%SNAKE2SWC convert a population of snakes to a single SWC file
xyz = vertcat(snks.vert)-1;                             % snk pnt coordinate
N = size(xyz,1);                                        % total number of snk pnt
swc = [(1:N)', 2*ones(N,1), xyz, ones(N,2)];            % neuron in SWC format
st = cumsum([0 arrayfun(@(x)size(x.vert,1), snks)]);	% snk's starting pnt number
for i = 1:length(snks)
    if snks(i).collide(1) == 0
        swc(st(i)+1,7) = -1;
    else
        prev = snks(i).collide;
        d = sum(bsxfun(@minus, snks(i).vert(1,:), snks(prev).vert).^2,2);
        [~,pos] = min(d);
        swc(st(i)+1,7) = st(prev)+pos;
    end
    n = size(snks(i).vert,1);
    swc(st(i)+2:st(i)+n,7) = st(i)+1:st(i)+n-1;
end
end

