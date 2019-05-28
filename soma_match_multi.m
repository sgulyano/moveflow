function soma1 = soma_match_multi( I1, I0, cens, rads, soma0, opt)
%SOMA_MATCH tracking soma location
% INPUT
%   I1 - current image
%   I0 - previous image
%   soma0 - previous soma
%   pdsoma2 - pair distance between soma in 2D at initial time step
%
if nargin < 8;  opt = struct();  end;
if ~isfield(opt,'filter_size');     opt.filter_size     = 20;   end;
if ~isfield(opt,'search_size');     opt.search_size     = 40;   end;
if ~isfield(opt,'iter_max');        opt.iter_max        = 100;  end;

filt_size = opt.filter_size;
search_size = opt.search_size;

n = size(soma0,1);
soma1 = zeros(n,3);

s = size(I1);
I0_pad = padarray(I0, [filt_size filt_size], 'both');
I1_pad = padarray(I1, [filt_size filt_size], 'both');

%% predict the new position using normalized cross correlation
cen_ncc = soma0(:,1:2);
score = zeros(1,n);
visited = zeros(1,n);
for i = 1:n
    cen0 = round(soma0(i,1:2));
    if all(cen0 == 0) || any(cen0 < -1 | cen0 > s([2 1]))
        visited(i) = -1;
        continue;
    end
    
    filter = im2double(I0_pad(cen0(2):cen0(2)+2*filt_size, cen0(1):cen0(1)+2*filt_size));

    topleft = max(cen0 - search_size, 1);
    bottomr = min(cen0 + search_size, s([2 1]));
    Icrop = im2double(I1_pad(topleft(2):bottomr(2)+2*filt_size, ...
            topleft(1):bottomr(1)+2*filt_size));

    C = normxcorr2(filter, Icrop);
    [ypeak, xpeak] = find(C==max(C(:)));
    yoffSet = ypeak - size(filter,1) + topleft(2);
    xoffSet = xpeak - size(filter,2) + topleft(1);

    if xoffSet < 1 || yoffSet < 1 || xoffSet > size(I1,2) || yoffSet > size(I1,1)
        continue
    end
    
    score(i) = max(C(:));
    cen_ncc(i,:) = [xoffSet yoffSet];
end

%% set new position using either NCC or candidates from Hough transform
[D,idx] = pdist2(cens, cen_ncc, 'euclidean', 'Smallest', n);
for i = 1:n^2
    [ii,jj] = find(D == min(D(:)), 1);
    if visited(jj)
        D(ii,jj) = inf;
        continue
    end
    old_rad = soma0(jj,3);

    if D(ii,jj) < old_rad
        soma1(jj,:) = [cens(idx(ii,jj),:) rads(idx(ii,jj))];
        D(idx==idx(ii,jj)) = inf;
    else
        try
        if score(jj) < 0.5 || sqrt(sum((soma0(jj,1:2) - cen_ncc(jj,:)).^2)) > 40
            soma1(jj,:) = zeros(1,3);
        else
            soma1(jj,:) = [cen_ncc(jj,:) old_rad];
        end
        catch
            keyboard;
        end
    end
    visited(jj) = idx(ii,jj);
end
end