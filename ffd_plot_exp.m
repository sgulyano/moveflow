function ffd_plot_exp( filename, V, swcs, track_time, soma_pos, dd_plane, opt )
%FFD_PLOT_BOTH draw results tracking in 3D over Depth using FFD and save to
%video

if nargin < 7;  opt = struct();  end
if ~isfield(opt,'is_gui');      opt.is_gui    = false;              end
if ~isfield(opt,'dx');          opt.dx        = 4;                  end
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end%0.112674;              end%231.282 / 1000;     end%
if ~isfield(opt,'maxdep');      opt.maxdep    = 27;                 end
if ~isfield(opt,'zflip');       opt.zflip     = [];                 end


configs = dd_plane.configs; 
spline = dd_plane.spline;


col1 = [80  200 80 ]/255; % ddaD color
col2 = [235 90  235]/255; % ddaE color

dx = opt.dx;

time = (0:size(soma_pos{1},1)-1) * opt.tpf;

sizeI = [size(V,1), size(V,2)];
Xmax = max(length(configs{track_time+1,1}), length(configs{track_time+1,2}))*dx;
Zmax = 40;%opt.maxdep*dx;

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
% get bounding box
[ii, jj] = find(I2D_mask);
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
I_mask1_crop = I_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
I_mask2_crop = I_mask2(yrange(1):yrange(2),xrange(1):xrange(2));


% init FFD
if isempty(opt.zflip)
    [Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx);
else
    [Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx, opt.zflip{track_time+1});
end
[sx, sy] = meshgrid((Ox_bs-1)*opt.pixelsize, (yrange-1)*opt.pixelsize);
sz = ones(2,1) * (Oz_bs+8) * opt.pixelsize;
    
% open video
if ~isempty(filename) && nargout ~= 1
    outputVideo = VideoWriter(filename, 'MPEG-4');
    outputVideo.Quality=95;
    outputVideo.FrameRate = 15;
    open(outputVideo);
end

if nargout ~= 1
    % create figure
    f = figure('color', 'w', 'units','pixels','position',[10 50 1024 576]);
%     f.InvertHardcopy = 'off';
    ha2 = axes('Units','normalized','Position',[0.65,0.05,0.29,0.9]);
%     ha2.XColor = 'w';
%     ha2.YColor = 'w';
%     ha2.ZColor = 'w';
    ha3 = axes('Units','normalized','Position',[0.05,0.05,0.55,0.9]);
%     ha3.XColor = 'w';
%     ha3.YColor = 'w';
%     ha3.ZColor = 'w';
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
    hffd_mod = surface(sx-(somaX-1)*opt.pixelsize, sy-(somaY-1)*opt.pixelsize, sz, 256+ones(2,1)*(1:size(sz,2)), ...
            'CDataMapping', 'direct', 'EdgeColor', 'interp','FaceAlpha',0.3);
    hffd_mod.AmbientStrength = 1.0;
%     hffd_mod.Visible = 'off';
    
    hold on;
    hmodel1 = patch(sf1);
    hmodel1.FaceColor = col1;
    hmodel1.EdgeColor = 'none';
    hmodel2 = patch(sf2);
    hmodel2.FaceColor = col2;
    hmodel2.EdgeColor = 'none';
%     patch(planebox([1 2 2 1]), planebox([3 3 4 4]), zeros(1,4), zeros(1,4),'FaceAlpha',.1);
    hold off;
    daspect([1 1 0.7])
    xlabel('X (\mum)');
    ylabel('Y (\mum)');
    zlabel('Z (\mum)');
    zlim([0 Zmax] * opt.pixelsize)
    % view([15 36]);
    view([15 30]);
    set(ha2,'Ydir','reverse');
    ha2.YDir = 'normal';
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
    hslice = surface(xx,yy,zz,255-V(:,:,track_time+1),...
            'FaceColor','texturemap',...
            'EdgeColor','none',...
            'CDataMapping','direct');
    hslice.AmbientStrength = 1.0;
    hslice.SpecularStrength = 0.0;
    colormap([gray(256);parula(size(sz,2))]);
    view(-35,45)
    set(ha3,'Ydir','reverse');
    ha3.YDir = 'normal';
    set(ha3,'Zdir','reverse');
    
    % draw FFD configuration
    hold on;
    hffd = surface(sx, sy, sz, 256+ones(2,1)*(1:size(sz,2)), 'CDataMapping', 'direct', 'EdgeColor', 'interp','FaceAlpha',0.1);
    hffd.AmbientStrength = 1.0;
%     hffd.Visible = 'off';
    
    sf1.vertices = bsxfun(@plus, sf1.vertices, [somaX-2, somaY-2, 0] * opt.pixelsize) ;
    sf2.vertices = bsxfun(@plus, sf2.vertices, [somaX-2, somaY-2, 0] * opt.pixelsize);
    hneu1 = patch(sf1);
    hneu1.FaceColor = col1;
    hneu1.EdgeColor = 'none';
    hneu2 = patch(sf2);
    hneu2.FaceColor = col2;
    hneu2.EdgeColor = 'none';
    
    hleft = scatter3(dd_plane.leftpos(2,1)*opt.pixelsize, ...
            dd_plane.leftpos(2,2)*opt.pixelsize, 0, 60, 'r*', 'LineWidth', 1.5);
    hright = scatter3(dd_plane.rightpos(2,1)*opt.pixelsize, ...
            dd_plane.rightpos(2,2)*opt.pixelsize, 0, 60, 'r*', 'LineWidth', 1.5);
%     hleft.Visible = 'off';
%     hright.Visible = 'off';
    
    hold off;
    box on
    axis equal
    view([0 90])
%     xlim([60 280])
%     ylim([20 size(V,1)*opt.pixelsize])
%     xlim([130 size(V,2)*opt.pixelsize])
%     ylim([7 70])
%     xlim([0 200])
%     ylim([10 82])
    zlim([0 Zmax*opt.pixelsize])
    camlight
    lighting gouraud
end
if ~isempty(filename) && nargout ~= 1
    [di, ~, ~] = fileparts(filename);
    saveas(f, fullfile(di, sprintf('%03d.png', track_time)));
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
        hneu1.Vertices = bsxfun(@minus, sf1.vertices, [2, 2, 0]) * opt.pixelsize;
        hneu1.Faces = sf1.faces;
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
        hneu2.Vertices = bsxfun(@minus, sf2.vertices, [2, 2, 0]) * opt.pixelsize;
        hneu2.Faces = sf2.faces;
        % relocate soma to (0,0,1)
        hmodel2.Vertices = bsxfun(@minus, sf2.vertices, [newsomaX, newsomaY, 0]) * opt.pixelsize;
        hmodel2.Faces = sf2.faces;
        
        % update image slice
        hslice.CData = 255-V(:,:,t);
    end
    
    Oz = (Oz - min(Oz));
%     disp(max(Oz)/dx);
    if nargout ~= 1
        % update FFD configuration
        hffd.XData = ones(2,1)*(Ox-1)*opt.pixelsize;
        hffd.YData = (newyrange'-1)*opt.pixelsize * ones(size(Ox));
        hffd.ZData = ones(2,1)*(Oz+8)*opt.pixelsize;
%         hffd.CData = ones(2,1)*(Oz+8) + 257;
        
        hffd_mod.XData = ones(2,1)*(Ox+1-newsomaX)*opt.pixelsize;
        hffd_mod.YData = (newyrange'+1-newsomaY)*opt.pixelsize * ones(size(Ox));
        hffd_mod.ZData = ones(2,1)*(Oz+8)*opt.pixelsize;
%         hffd_mod.CData = ones(2,1)*(Oz+8) + 257;
        pause(0.1);
    end
    
    hleft.XData = dd_plane.leftpos(t,1)*opt.pixelsize;
    hleft.YData = dd_plane.leftpos(t,2)*opt.pixelsize;
    hright.XData = dd_plane.rightpos(t,1)*opt.pixelsize;
    hright.YData = dd_plane.rightpos(t,2)*opt.pixelsize;
    
%     axes(ha3);
%     title(sprintf('Frame %d/%d', t-1, size(V,3)-1));
    title(sprintf('Time %.3f s, Frame %d/%d', time(t), t-1, size(V,3)-1));
    
    % save to video
    if ~isempty(filename) && nargout ~= 1
        
        saveas(f, fullfile(di, sprintf('%03d.png', t-1)));
        writeVideo(outputVideo, getframe(gcf));
    end
end
if ~isempty(filename) && nargout ~= 1
    close(outputVideo);
end
pause;
delete(f);
end

