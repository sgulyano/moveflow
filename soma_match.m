function [centroid, rad, flag] = soma_match( I, Iprev, centers, radii, old_centroid, rad, opt)
%SOMA_MATCH tracking soma location using normalized cross-correlation

if nargin < 7;  opt = struct();  end;
if ~isfield(opt,'filter_size');     opt.filter_size     = 20;   end;
if ~isfield(opt,'search_size');     opt.search_size     = 50;   end;

s = size(I);
if all(old_centroid == 0) || any(old_centroid < -1 | old_centroid > s([2 1]))
    centroid = [0 0];
    flag = true;
    return;
end

old_centroid = round(old_centroid);

filter_size = opt.filter_size;
search_size = opt.search_size;

Iprev_pad = padarray(Iprev, [filter_size filter_size], 'both');
filter = im2double(Iprev_pad(old_centroid(2):old_centroid(2)+2*filter_size, ...
        old_centroid(1):old_centroid(1)+2*filter_size));

I_pad = padarray(I, [filter_size filter_size], 'both');


topleft = max(old_centroid - search_size, 1);
bottomr = min(old_centroid + search_size, s([2 1]));

Icrop = im2double(I_pad(topleft(2):bottomr(2)+2*filter_size, ...
        topleft(1):bottomr(1)+2*filter_size));

C = normxcorr2(filter, Icrop);
valpeak = max(C(:));
[ypeak, xpeak] = find(C==max(C(:)));
yoffSet = ypeak-size(filter,1)+topleft(2);
xoffSet = xpeak-size(filter,2)+topleft(1);

if xoffSet < 1 || yoffSet < 1 || xoffSet > size(I,2) || yoffSet > size(I,1)
    centroid = [xoffSet yoffSet];
    flag = true;
    return;
end

% figure; imshow(I);
% hold on;
% scatter(xoffSet, yoffSet, 'rx');

[ds,pos] = min(sqrt(sum(bsxfun(@minus, centers, [xoffSet, yoffSet]).^2,2)));
if ds < rad
    centroid = centers(pos,:);
    rad = radii(pos);
else
    centroid = [xoffSet yoffSet];
end

% if valpeak < 0.6 || sqrt(sum((old_centroid-centroid).^2)) > 40%ds > 40,
%     flag = true;
%     keyboard;
%     return;
% else
    flag = false;
% end

% figure(2); imshow(I);
% viscircles(centers, radii,'EdgeColor','r');

end

