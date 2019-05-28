function [Ox, Oz] = plane_init_from_config( config, somaX, dx )
%FFD_INIT_FROM_CONFIG given configuration of dendrite fold as list of
%orientation and control point distribution wrt soma and X-coordinate of
%soma with spacing dx, return the FFD control line/grid

d1 = cumsum(config{1});
d2 = cumsum(config{2});

% if nargin < 4
%     zflip = ones(1, length(d1) + length(d2) + 1);
% end
% zflip_cell = cell(1,2);
% zflip_cell{1} = sign(zflip(length(d1):-1:1));
% zflip_cell{2} = sign(zflip(length(d1)+2:end));


Ox = [fliplr(-d1) 0 d2] + somaX;
if nargout > 1
    Oz = [fliplr(cumsum(dx * sqrt(1 - (config{1} / dx).^2))), 0, cumsum(dx * sqrt(1 - (config{2} / dx).^2))];
end
end

