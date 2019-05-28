function [ Ox, Oz ] = ffd_init_from_config( config, somaX, dx, zflip )
%FFD_INIT_FROM_CONFIG given configuration of dendrite fold as list of
%orientation and control point distribution wrt soma and X-coordinate of
%soma with spacing dx, return the FFD control line/grid

ang1 = cumsum(config{1});
ang2 = cumsum(config{2});

Ox = [fliplr(-cumsum(cos(ang1)*dx)) 0 cumsum(cos(ang2)*dx)] + somaX;
if nargout > 1
    if nargin < 4
        zflip = ones(1, length(ang1) + length(ang2) + 1);
    end
    if iscell(zflip)
        zflip_cell = zflip;
    else
        zflip_cell = cell(1,2);
        zflip_cell{1} = sign(zflip(length(ang1):-1:1));
        zflip_cell{2} = sign(zflip(length(ang1)+2:end));
    end
    
    Oz = [fliplr(cumsum( abs(sin(ang1)*dx).*zflip_cell{1}  )) 0 cumsum( abs(sin(ang2)*dx).*zflip_cell{2} )]/2;
%     Oz2 = [fliplr(cumsum(abs(sin(ang1)*dx))) 0 cumsum(abs(sin(ang2)*dx))];
end
end

