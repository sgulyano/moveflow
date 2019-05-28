function [ cutindex ] = check_intersection( swcdata, cutline )
%CHECK_INTERSECTION find edges in swcdata that intersects with cutline

% let endpoint of cutline be P and Q and let endpoint of dendrite edges be
% U and V. Check intersection by sign(PQxQU) ~= sign(PQxQV) and 
% sign(UVxVP) ~= sign(UVxVQ)

U = swcdata(:,3:4);
idx = swcdata(:,7)>0;
prev = swcdata(idx,7);
V = swcdata(:,3:4); V(idx,:) = swcdata(prev,3:4);

PQ = cutline(2,:) - cutline(1,:);
QU = bsxfun(@minus, U, cutline(2,:));
QV = bsxfun(@minus, V, cutline(2,:));

PQxQU = PQ(1)*QU(:,2) - PQ(2)*QU(:,1);
PQxQV = PQ(1)*QV(:,2) - PQ(2)*QV(:,1);

UV = V - U;
VP = bsxfun(@minus, cutline(1,:), V);
VQ = -QV;
UVxVP = UV(:,1).*VP(:,2) - UV(:,2).*VP(:,1);
UVxVQ = UV(:,1).*VQ(:,2) - UV(:,2).*VQ(:,1);

cutindex = (sign(PQxQU) ~= sign(PQxQV)) & (sign(UVxVP) ~= sign(UVxVQ));


end

