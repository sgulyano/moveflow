function soma1 = soma_match_multi( I1, I0, cens, rads, metrics, soma0, pdsoma2, opt)
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

%% First, find candidate soma location from:
cand_soma = cell(1,size(soma0,1));

% a) NCC-predicted position
s = size(I1);
I0_pad = padarray(I0, [filt_size filt_size], 'both');
I1_pad = padarray(I1, [filt_size filt_size], 'both');

cen_ncc = zeros(n, 2);
rad_ncc = zeros(n, 1);
for i = 1:n
    cen0 = round(soma0(i,1:2));
    if all(cen0 == 0) || any(cen0 < -1 | cen0 > s([2 1]))
        continue;
    end
    
    filter = im2double(I0_pad(cen0(2):cen0(2)+2*filt_size, cen0(1):cen0(1)+2*filt_size));

    topleft = max(cen0 - search_size, 1);
    bottomr = min(cen0 + search_size, s([2 1]));
    I1_crop = im2double(I1_pad(topleft(2):bottomr(2)+2*filt_size, ...
            topleft(1):bottomr(1)+2*filt_size));

    C = normxcorr2(filter, I1_crop);
    [ypeak, xpeak] = find(C==max(C(:)));
    yoffSet = ypeak - size(filter,1) + topleft(2);
    xoffSet = xpeak - size(filter,2) + topleft(1);

    if xoffSet < 1 || yoffSet < 1 || xoffSet > size(I1,2) || yoffSet > size(I1,1)
        continue
    end
    
    cen_ncc(i,:) = [xoffSet yoffSet];
    rad_ncc(i) = soma0(i,3);
    cand_soma{i} = size(cens,1);
end


% b) Previous position
cens = [cens; cen_ncc];
rads = [rads; rad_ncc];
metrics = [metrics; zeros(n,1)];

% c) Circles (from Hough) close to previous or NCC-predicted positions
for i = 1:n
    cand_soma{i} = unique([find(sum(bsxfun(@minus, cens, cen_ncc(i,:)).^2,2) < opt.search_size.^2); ...
            find(sum(bsxfun(@minus, cens, soma0(i,1:2)).^2,2) < opt.search_size.^2)]);
end

%% Second, optimize obj func using ICM, where candidates are states
% initial state is the previous position
state = cellfun(@(x)x(end), cand_soma);
Emin = get_energy(I1, cens, rads, metrics, pdsoma2, state);
% ICM starts
for iter = 1:opt.iter_max
    % pick random variable
    ii = randi(n);
    state_tmp = state;
    % find the optimal state while fixed other variables
    for cand = cand_soma{ii}'
        state_tmp(ii) = cand;
        E = get_energy(I1, cens, rads, metrics, pdsoma2, state_tmp);
        if E < Emin
            Emin = E;
            state = state_tmp;
        end
    end
end

soma1 = [cens(state,:), rads(state)];

end

function E = get_energy(I, cens, rads, metrics, pdsoma2, state)
    cen = cens(state,:);
    rad = rads(state);
    ds2 = 0;
    for i = 1:2
        ds2 = ds2 + sum(sum(abs(tril(bsxfun(@minus, cen(:,i), cen(:,i)') - pdsoma2(:,:,i)))));
    end
    
    [xx, yy] = meshgrid(1:size(I,2),1:size(I,1));
    mask = false(size(I));
    for i = 1:size(cen,1)
        mask(sqrt((xx-cen(i,1)).^2 + (yy-cen(i,2)).^2) <= rad(i)) = true;
    end    
    
    keyboard;
    
    E = exp(-ds2/20) + 10*sum(I(mask))/(255*sqrt(numel(mask)));
end