function [Vseed, seedp] = IV_findseed( Iout, mask, varargin )
%FINDSEEDPOINT find seed volume and seed points in XYZ-coordinate from
%preprocess volume using local maxima of OOF response. Hessian based
%technique like Frangi give false positive because of clutter neurite.

defaultoptions = struct('gap', 3);
options = get_options(defaultoptions, varargin);

REPORT = false;
rect = [30, 30, 200, 200];

%% Find candidate points
if ismatrix(Iout)
    [fx,fy] = AM_gradient(Iout);
else
    [fx,fy,~] = AM_gradient(Iout);
end
if nnz(mask) > 1
    seeds = Skeleton3D(mask);
else
    OOF_thr = mask;
    seedsx = fx(:,1:end-1,:)>0 & fx(:,2:end,:)<0 & Iout(:,1:end-1,:)>OOF_thr;
    seedsy = fy(1:end-1,:,:)>0 & fy(2:end,:,:)<0 & Iout(1:end-1,:,:)>OOF_thr;
    seeds = seedsx(1:end-1,:,:) | seedsy(:,1:end-1,:);
    seeds = padarray(seeds,[1 1],'symmetric','post');
end


if REPORT
    tmp = imcrop(max(seeds, [], 3), rect);
    figure(8), imshow(tmp, [], 'Border','tight');
end

%% keep local maxima of OOF response as seed point
gap = -options.gap:options.gap;
[x,y,z] = meshgrid(gap,gap,gap);
% [x,y,z] = meshgrid(-3:3,-3:3,-3:3);
x = x(:); y = y(:); z = z(:);
if ismatrix(Iout)
    s = [size(Iout), 1];
else
    s = size(Iout);
end
max_size = ones(numel(x),1) * s([2,1,3]);
p_ind = find(seeds > 0);
[py, px, pz] = ind2sub(s, p_ind);
keep = false(1,length(p_ind));
parfor i = 1:length(p_ind)
    % check each candidate against other candidates in 27-neighborhood
    nb = [px(i)+x, py(i)+y, pz(i)+z];
    cond = all(nb>=1,2) & all(nb<=max_size,2);
    idx = sub2ind(s, py(i)+y(cond), px(i)+x(cond), pz(i)+z(cond));
    idx = idx(seeds(idx));      % max among local seed candidates
    [~,maxpos] = max(Iout(idx));
    % keep the local maximum of OOF response
    if p_ind(i) == idx(maxpos)
        keep(i) = true;
    end
end

% get XYZ-coordinate of seed point and seed volume
Vseed = false(s);
Vseed(p_ind(keep)) = true;
seedp = [px(keep), py(keep), pz(keep)];

if REPORT
    tmp = imcrop(max(Vseed, [], 3), rect);
    figure(9), imshow(tmp, [], 'Border','tight');
%     keyboard;
end
end