function ffd_plot( filename, V, swc, track_time, soma_pos, configs, spline, dx, opt )
%FFD_PLOT draw results tracking in 3D over Depth using FFD and save to
%video

if nargin < 9;  opt = struct();  end;
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end;
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end;

if nargin < 8
    dx = 4;
end

time = (0:size(soma_pos,1)-1) * opt.tpf;

sizeI = [size(V,1), size(V,2)];
Zmax = max(length(configs{track_time+1,1}), length(configs{track_time+1,2}))*dx;
Zmin = 0;%length(configs{1,2})*dx;

% create baseline mask
I_mask = swc2pixel( swc, sizeI );
I2D_mask = imdilate(I_mask, ones(3));
% get bounding box
[ii, jj] = find(I2D_mask);
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
I_mask_crop = I_mask(yrange(1):yrange(2),xrange(1):xrange(2));
I2D_mask_crop = I2D_mask(yrange(1):yrange(2),xrange(1):xrange(2));
% locate soma
somaX = swc(1,3);
somaY = swc(1,4);

% init FFD
[Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx);

% open video
if ~isempty(filename)
    outputVideo = VideoWriter(filename, 'MPEG-4');
    open(outputVideo);
end

% create figure
f = figure('Visible','off','Position',[160,200,1050,785]);
ha1 = axes('Units','normalized','Position',[0.05,0.55,0.6,0.4]);
ha2 = axes('Units','normalized','Position',[0.7,0.55,0.25,0.4]);
ha3 = axes('Units','normalized','Position',[0.05,0.05,0.9,0.45]);

% draw image and highlight dendrite
axes(ha1);
himg = imshow(repmat(V(:,:,track_time+1),1,1,3));
green = cat(3, zeros(sizeI), ones(sizeI), zeros(sizeI));
hold on;
hmask = imshow(green);
set(hmask, 'AlphaData', I2D_mask*.2);
hold off;

% transform dendrite
Tx = ffd_interpolate(Ox_bs, spline);
Tz = ffd_interpolate(Oz_bs, spline)+1;
[cx, cy] = meshgrid(round(Tx), yrange(1):yrange(2));
[cz, ~] = meshgrid(round(Tz), yrange(1):yrange(2));
cx = cx(:); cy = cy(:); cz = cz(:);

% create 3D surface of dendrite morphology
sizeI3D = [sizeI 1+range(cz)];
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(I_mask_crop), cx(I_mask_crop), cz(I_mask_crop))) = 1;
I3D = padarray(I3D, [2 2 2], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
sf = isosurface(I3D, .5);
% relocate soma to (0,0,1)
sf.vertices = bsxfun(@minus, sf.vertices, [somaX, somaY, 3-min(cz)]) * opt.pixelsize;

% draw 3D model of dendrite
axes(ha2);cla;
hmodel = patch(sf);
hmodel.FaceColor = 'red';
hmodel.EdgeColor = 'none';
daspect([1 1 1])
xlabel('X (\mum)');
ylabel('Y (\mum)');
zlabel('Z (\mum)');
planebox = [-Zmax Zmax yrange-somaY+[-10 10] Zmin-2 Zmax] * opt.pixelsize;
axis([-Zmax Zmax yrange-somaY+[-10 10] Zmin-2 Zmax] * opt.pixelsize)
zlim([Zmin-2 Zmax] * opt.pixelsize)
view([15 36]);
set(ha2,'Ydir','reverse');
set(ha2,'Zdir','reverse');
hold on;
patch(planebox([1 2 2 1]), planebox([3 3 4 4]), -2*opt.pixelsize*ones(1,4), zeros(1,4),'FaceAlpha',.1);
hold off;

box on
camlight 
lighting gouraud

% draw image slice
[xx, yy, zz] = meshgrid((0:sizeI(2)-1)*opt.pixelsize, (0:sizeI(1)-1)*opt.pixelsize, Zmax*opt.pixelsize);
axes(ha3);
title(sprintf('Time %.3f s, Frame %d', time(track_time+1), track_time));
xlabel('X (\mum)');
ylabel('Y (\mum)');
zlabel('Z (\mum)');
hslice = surface(xx,yy,zz,V(:,:,track_time+1),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
colormap([gray(256);parula(Zmin+Zmax)]);
view(-35,45)
set(ha3,'Ydir','reverse');
set(ha3,'Zdir','reverse');
% draw FFD configuration
hold on;
[sx, sy] = meshgrid((Ox_bs-1)*opt.pixelsize, (yrange-1)*opt.pixelsize);
sz = ones(2,1) * Oz_bs*opt.pixelsize;
hffd = surface(sx, sy, sz, sz+257+Zmin, 'CDataMapping', 'direct', 'EdgeColor', 'interp','FaceAlpha',0.5);
hold off;
axis equal
box on
view([-19 25])
zlim([0 Zmax*opt.pixelsize])

% save to video
f.Visible = 'on';
if ~isempty(filename)
    writeVideo(outputVideo, getframe(f));
end

%% create animation
for t = track_time+2:size(V,3)
    if isempty(configs{t,1})
        break;
    end
    [Ox,Oz] = ffd_init_from_config(configs(t,:), soma_pos(t,1), dx);
    % transform dendrite in 2D
    Tx = ffd_interpolate(Ox, spline);
    newyrange = yrange - round(somaY - soma_pos(t,2));
    [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
    cx = cx(:); cy = cy(:);
    I_mask_t = zeros(sizeI);
    idxIn = cy >= 1 & cy <= sizeI(1) & cx >= 1 & cx <= sizeI(2);
    I_mask_t(sub2ind(sizeI, cy(I2D_mask_crop(:)&idxIn), cx(I2D_mask_crop(:)&idxIn))) = 1;
    % update image and dendrite highlight
    set(himg, 'CData', repmat(V(:,:,t),1,1,3));
    title(sprintf('Time %.3f s, Frame %d', time(t), t-1));
    set(hmask, 'AlphaData', I_mask_t*.2);
    
    % transform dendrite in Z-axis
    Tz = ffd_interpolate(Oz, spline)+1;
    [cz, ~] = meshgrid(round(Tz) - min(round(Tz))+1, newyrange(1):newyrange(2));
    cz = cz(:);
    % update 3D dendrite morphology model
    sizeI3D = [sizeI 1+range(cz)];
    I3D = zeros(sizeI3D);
    I3D(sub2ind(sizeI3D, cy(I_mask_crop(:)&idxIn), cx(I_mask_crop(:)&idxIn), cz(I_mask_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 2], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    sf = isosurface(I3D, .5);
    % relocate soma to (0,0,1)
    hmodel.Vertices = bsxfun(@minus, sf.vertices, [soma_pos(t,1:2), 3-min(cz)]) * opt.pixelsize;
    hmodel.Faces = sf.faces;

    % update image slice
    hslice.CData = V(:,:,t);
    % update FFD configuration
    hffd.XData = ones(2,1)*(Ox-1)*opt.pixelsize;
    hffd.YData = (newyrange'-1)*opt.pixelsize * ones(size(Ox));
    hffd.ZData = ones(2,1)*Oz*opt.pixelsize;
    hffd.CData = ones(2,1)*Oz + 257 + Zmin;
    
    % save to video
    pause(0.1);
    if ~isempty(filename)
        writeVideo(outputVideo, getframe(gcf));
    end
end
if ~isempty(filename)
    close(outputVideo);
end
pause;
delete(f);
end

