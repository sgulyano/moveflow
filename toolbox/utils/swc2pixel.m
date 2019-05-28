function img = swc2pixel( swc, Isize )
%SWC2PIXEL draw the trace in swc format in pixel
if length(Isize)==2, Isize = [Isize inf]; end
swc(:,3:5) = round(swc(:,3:5));
X = cell(1,size(swc,1));
n = size(swc,1);
% for every line, find pixels that form the line.
idxOut = ~all(swc(:,3:5) >= 1 & swc(:,3:5) <= ones(n,1)*Isize([2 1 3]), 2);
for i = 1:n
    prev = swc(i,7);
    if prev < 1 || idxOut(i) || idxOut(prev), continue; end
    % find valid line (ignore lines with root node, out-of-bound endpoint)
    X{i} = getpixelline(swc(i,3:4), swc(prev,3:4));
end
P = round(vertcat(X{:}));
% convert to binary image
img = false(Isize(1:2));
if isempty(P)
    return;
end
idx = sub2ind(Isize(1:2), P(:,2), P(:,1));
img(idx) = true;

    function X = getpixelline(P1, P2)
    % find pixels that form a line bewteen P1 and P2
        dP = P2-P1;
        len = ceil(sqrt(sum(dP.^2)));
        steps = 0:1/len:1;
        X = ones(length(steps),1)*P1 + steps'*dP;
    end
end
