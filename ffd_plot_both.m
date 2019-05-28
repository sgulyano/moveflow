function masks = ffd_plot_both( filename, V, swcs, track_time, soma_pos, configs, spline, opt )
%FFD_PLOT_BOTH draw results tracking in 3D over Depth using FFD and save to
%video

if nargin < 8;  opt = struct();  end;
if ~isfield(opt,'is_gui');      opt.is_gui    = false;              end;
if ~isfield(opt,'dx');          opt.dx        = 4;                  end;
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end;
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end;
if ~isfield(opt,'maxdep');      opt.maxdep    = 27;                 end;
if ~isfield(opt,'zflip');       opt.zflip     = [];                 end;

if nargout == 1,
    masks = cell(2, size(V,3));
end

col1 = [235 90  235]/255; % ddaD color
col2 = [80  200 80 ]/255; % ddaE color

dx = opt.dx;

time = (0:size(soma_pos{1},1)-1) * opt.tpf;

sizeI = [size(V,1), size(V,2)];
Xmax = max(length(configs{track_time+1,1}), length(configs{track_time+1,2}))*dx;
Zmax = opt.maxdep*dx;

% locate soma
ddad_somaX = round(swcs{1}(1,3));
ddad_somaY = round(swcs{1}(1,4));
ddae_somaX = round(swcs{2}(1,3));
ddae_somaY = round(swcs{2}(1,4));
somaX = round((ddad_somaX + ddae_somaX)/2);
somaY = round((ddad_somaY + ddae_somaY)/2);

% create baseline mask
I_mask1 = swc2pixel( swcs{1}, sizeI );
I_mask2 = swc2pixel( swcs{2}, sizeI );
I2D_mask = imdilate(I_mask1 | I_mask2, ones(3));
I2D_mask1 = imdilate(I_mask1, ones(3));
I2D_mask2 = imdilate(I_mask2, ones(3));

[xx, yy] = meshgrid(1:size(I2D_mask,2), 1:size(I2D_mask,1));
I2D_mask1((xx - ddad_somaX).^2 + (yy - ddad_somaY).^2 < swcs{1}(1,6)^2) = 1;
I2D_mask2((xx - ddae_somaX).^2 + (yy - ddae_somaY).^2 < swcs{2}(1,6)^2) = 1;
% get bounding box
[ii, jj] = find(I2D_mask);
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
I_mask1_crop = I_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
I_mask2_crop = I_mask2(yrange(1):yrange(2),xrange(1):xrange(2));
I2D_mask1_crop = I2D_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
I2D_mask2_crop = I2D_mask2(yrange(1):yrange(2),xrange(1):xrange(2));

if nargout == 1,
    masks{1, track_time+1} = I2D_mask1;
    masks{2, track_time+1} = I2D_mask2;
end

% init FFD
if isempty(opt.zflip)
    [Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx);
else
    [Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx, opt.zflip{track_time+1});
end

% open video
if ~isempty(filename) && nargout ~= 1
    outputVideo = VideoWriter(filename, 'MPEG-4');
    open(outputVideo);
end

if nargout ~= 1
    % create figure
    f = figure('Visible','off','Position',[160,200,1050,785]);
    ha1 = axes('Units','normalized','Position',[0.05,0.55,0.6,0.4]);
    ha2 = axes('Units','normalized','Position',[0.7,0.55,0.25,0.4]);
    ha3 = axes('Units','normalized','Position',[0.05,0.05,0.9,0.45]);

    % draw image and highlight dendrite
    axes(ha1);
    himg = imshow(repmat(V(:,:,track_time+1),1,1,3));
    hold on;
    hmask1 = imshow(repmat(reshape(col1, [1 1 3]), sizeI));
    set(hmask1, 'AlphaData', I2D_mask1*.2);
    hmask2 = imshow(repmat(reshape(col2, [1 1 3]), sizeI));
    set(hmask2, 'AlphaData', I2D_mask2*.2);
    hold off;
else
    if opt.is_gui
        fwait = waitbar(0,'Computing Mask of Dendrites');
    end
end

% transform dendrite
Tx = ffd_interpolate(Ox_bs, spline);
Tz = ffd_interpolate(Oz_bs, spline)+1;
[cx, cy] = meshgrid(round(Tx), yrange(1):yrange(2));
[cz, ~] = meshgrid(round(Tz), yrange(1):yrange(2));
ddad_somaZ = cz(ddad_somaY - yrange(1) + 1, ddad_somaX - xrange(1) + 1);
ddae_somaZ = cz(ddae_somaY - yrange(1) + 1, ddae_somaX - xrange(1) + 1);
cx = cx(:); cy = cy(:); cz = cz(:);

% create 3D surface of ddaD
sizeI3D = [sizeI 1+range(cz)];
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(I_mask1_crop), cx(I_mask1_crop), cz(I_mask1_crop))) = 1;
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
[xx3, yy3, zz3] = meshgrid(1:size(I3D,2), 1:size(I3D,1), 1:size(I3D,3));
I3D((xx3 - ddad_somaX - 2).^2 + (yy3 - ddad_somaY - 2).^2 + (zz3 - ddad_somaZ - 7).^2 < swcs{1}(1,6)^2) = 1;
sf1 = isosurface(I3D, .5);
% relocate soma to (0,0,1)
sf1.vertices = bsxfun(@minus, sf1.vertices, [somaX, somaY, 0]) * opt.pixelsize;

% create 3D surface of ddaE
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(I_mask2_crop), cx(I_mask2_crop), cz(I_mask2_crop))) = 1;
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
I3D((xx3 - ddae_somaX - 2).^2 + (yy3 - ddae_somaY - 2).^2 + (zz3 - ddae_somaZ - 7).^2 < swcs{2}(1,6)^2) = 1;
sf2 = isosurface(I3D, .5);
% relocate soma to (0,0,1)
sf2.vertices = bsxfun(@minus, sf2.vertices, [somaX, somaY, 0]) * opt.pixelsize;

% draw 3D model of dendrite
planebox = [-Xmax Xmax yrange-somaY+[-10 10] 0 Zmax] * opt.pixelsize;

if nargout ~= 1
    axes(ha2);cla; axis(planebox);
    hmodel1 = patch(sf1);
    hmodel1.FaceColor = col1;
    hmodel1.EdgeColor = 'none';
    hold on;
    hmodel2 = patch(sf2);
    hmodel2.FaceColor = col2;
    hmodel2.EdgeColor = 'none';
    patch(planebox([1 2 2 1]), planebox([3 3 4 4]), zeros(1,4), zeros(1,4),'FaceAlpha',.1);
    hold off;
    daspect([1 1 1])
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    zlabel('Z (\mum)');
    zlim([0 Zmax] * opt.pixelsize)
    % view([15 36]);
    view([0 0]);
    set(ha2,'Ydir','reverse');
    set(ha2,'Zdir','reverse');
    box on
    camlight 
    lighting gouraud

    % draw image slice
    [xx, yy, zz] = meshgrid((0:sizeI(2)-1)*opt.pixelsize, (0:sizeI(1)-1)*opt.pixelsize, Zmax * opt.pixelsize);
    axes(ha3);
    title(sprintf('Time %.3f s, Frame %d/%d', time(track_time+1), track_time, size(V,3)-1));
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    zlabel('Z (\mum)');
    hslice = surface(xx,yy,zz,V(:,:,track_time+1),...
            'FaceColor','texturemap',...
            'EdgeColor','none',...
            'CDataMapping','direct');
    colormap([gray(256);parula(Zmax)]);
    view(-35,45)
    set(ha3,'Ydir','reverse');
    set(ha3,'Zdir','reverse');
    % draw FFD configuration
    hold on;
    [sx, sy] = meshgrid((Ox_bs-1)*opt.pixelsize, (yrange-1)*opt.pixelsize);
    sz = ones(2,1) * Oz_bs*opt.pixelsize;
    hffd = surface(sx, sy, sz, sz+257, 'CDataMapping', 'direct', 'EdgeColor', 'interp','FaceAlpha',0.5);
    hold off;
    box on
    axis equal
    view([-19 25])
    zlim([0 Zmax*opt.pixelsize])

    % save to video
    f.Visible = 'on';
end
if ~isempty(filename) && nargout ~= 1
    writeVideo(outputVideo, getframe(f));
end


%% create animation
for t = track_time+2:size(V,3)
    if isempty(configs{t,1})
        break;
    end
    ddadX = soma_pos{1}(t,1);
    ddadY = soma_pos{1}(t,2);
    ddaeX = soma_pos{2}(t,1);
    ddaeY = soma_pos{2}(t,2);
    newsomaX = (ddadX + ddaeX)/2;
    newsomaY = (ddadY + ddaeY)/2;
    dsomaY = round(somaY - (ddadY+ddaeY)/2);
    if isempty(opt.zflip)
        [Ox,Oz] = ffd_init_from_config(configs(t,:), newsomaX, dx);
    else
        [Ox,Oz] = ffd_init_from_config(configs(t,:), newsomaX, dx, opt.zflip{t});
    end
    % transform dendrite in 2D
    Tx = ffd_interpolate(Ox, spline);
    newyrange = yrange - dsomaY;
    [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
    cx = cx(:); cy = cy(:);
    idxIn = cy >= 1 & cy <= sizeI(1) & cx >= 1 & cx <= sizeI(2);
    
    I_mask_t1 = zeros(sizeI);
    I_mask_t1(sub2ind(sizeI, cy(I2D_mask1_crop(:)&idxIn), cx(I2D_mask1_crop(:)&idxIn))) = 1;
    I_mask_t2 = zeros(sizeI);
    I_mask_t2(sub2ind(sizeI, cy(I2D_mask2_crop(:)&idxIn), cx(I2D_mask2_crop(:)&idxIn))) = 1;
    if nargout ~= 1
        % update image and dendrite highlight
        set(himg, 'CData', repmat(V(:,:,t),1,1,3));
        title(sprintf('Time %.3f s, Frame %d/%d', time(t), t-1, size(V,3)-1));
        set(hmask1, 'AlphaData', I_mask_t1*.2);
        set(hmask2, 'AlphaData', I_mask_t2*.2);
    else
        masks{1, t} = I_mask_t1;
        masks{2, t} = I_mask_t2;
        if opt.is_gui
            waitbar((t-track_time-1)/(size(V,3)-track_time-1),fwait,'Computing Mask of Dendrites');
        end
    end
    
    % transform dendrite in Z-axis
    Tz = ffd_interpolate(Oz, spline)+1;
    [cz, ~] = meshgrid(round(Tz) - min(round(Tz))+1, newyrange(1):newyrange(2));
    ddadZ = cz(ddad_somaY - yrange(1) + 1, ddad_somaX - xrange(1) + 1);
    ddaeZ = cz(ddae_somaY - yrange(1) + 1, ddae_somaX - xrange(1) + 1);
    cz = cz(:);
    sizeI3D = [sizeI 1+range(cz)];
    % update 3D dendrite ddaD
    I3D = zeros(sizeI3D);
    I3D(sub2ind(sizeI3D, cy(I_mask1_crop(:)&idxIn), cx(I_mask1_crop(:)&idxIn), cz(I_mask1_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    [xx3, yy3, zz3] = meshgrid(1:size(I3D,2), 1:size(I3D,1), 1:size(I3D,3));
    I3D((xx3 - ddadX - 2).^2 + (yy3 - ddadY - 2).^2 + (zz3 - ddadZ - 7).^2 < swcs{1}(1,6)^2) = 1;
    sf1 = isosurface(I3D, .5);
    if nargout ~= 1
        % relocate soma to (0,0,1)
        hmodel1.Vertices = bsxfun(@minus, sf1.vertices, [newsomaX, newsomaY, 0]) * opt.pixelsize;
        hmodel1.Faces = sf1.faces;
    end
    % update 3D dendrite ddaE
    I3D = zeros(sizeI3D);
    I3D(sub2ind(sizeI3D, cy(I_mask2_crop(:)&idxIn), cx(I_mask2_crop(:)&idxIn), cz(I_mask2_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    I3D((xx3 - ddaeX - 2).^2 + (yy3 - ddaeY - 2).^2 + (zz3 - ddaeZ - 7).^2 < swcs{2}(1,6)^2) = 1;    
    sf2 = isosurface(I3D, .5);
    if nargout ~= 1
        % relocate soma to (0,0,1)
        hmodel2.Vertices = bsxfun(@minus, sf2.vertices, [newsomaX, newsomaY, 0]) * opt.pixelsize;
        hmodel2.Faces = sf2.faces;
        
        % update image slice
        hslice.CData = V(:,:,t);
    end
    
    Oz = (Oz - min(Oz));
%     disp(max(Oz)/dx);
    if nargout ~= 1
        % update FFD configuration
        hffd.XData = ones(2,1)*(Ox-1)*opt.pixelsize;
        hffd.YData = (newyrange'-1)*opt.pixelsize * ones(size(Ox));
        hffd.ZData = ones(2,1)*Oz*opt.pixelsize;
        hffd.CData = ones(2,1)*Oz + 257;
        pause(0.1);
    end
    
    % save to video
    if ~isempty(filename) && nargout ~= 1
        writeVideo(outputVideo, getframe(gcf));
    end
end
if ~isempty(filename) && nargout ~= 1
    close(outputVideo);
end
if nargout < 1,
    if opt.is_gui
        uiwait(msgbox('Finish'));
    else
        pause;
    end
    delete(f);
else
    if opt.is_gui
        delete(fwait);
    end
end
end

