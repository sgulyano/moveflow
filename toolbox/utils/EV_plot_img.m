function hs1 = EV_plot_img( im, data, col, varargin )
%EV_PLOT_IMG plot image with overlaid trace

if nargin < 3
    col = 'c';
end

if nargout == 0 && ~isempty(im),
    imshow(im);
end

if isempty(data)
    hs1 = [];
    return;
end
nb = data(:,[1 7]);
nb(any(nb < 1,2),:) = [];

G = sparse(nb(:,1), nb(:,2), 1, size(data,1), size(data,1));
G = G + G';
[ccnum, cc] = graphconncomp(G);

if ccnum == 1
    xx = [data(nb(:,1),3), data(nb(:,2),3), nan(size(nb,1),1)]';
    yy = [data(nb(:,1),4), data(nb(:,2),4), nan(size(nb,1),1)]';

    hold on;
    hs1 = plot(xx(:), yy(:), col, varargin{:});
    hold off;
else
    col = jet(max(ccnum,5));
    %col = jet(max(20, length(sizes)+6));
    %col = col(6:end,:);
    hs1 = [];
    hold on;
    for j = 1:ccnum
        idx = all(cc(nb)==j,2);
        xx = [data(nb(idx,1),3), data(nb(idx,2),3), nan(sum(idx),1)]';
        yy = [data(nb(idx,1),4), data(nb(idx,2),4), nan(sum(idx),1)]';
        hs1(j) = plot(xx(:), yy(:), 'LineWidth', 1, 'Color', col(j,:));
    end
    hold off;
end
end

