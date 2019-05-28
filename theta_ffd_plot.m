function theta_ffd_plot( filename, V, swc, track_time, soma_pos, soma_adj, configs, spline, opt )
%FFD_PLOT draw results tracking in 3D over Depth using FFD and save to
%video

if nargin < 9;  opt = struct();  end;
if ~isfield(opt,'dx');          opt.dx = 4;                         end;
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end;
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end;
if ~isfield(opt,'tsem');        opt.tsem   = [];                    end;
if ~isfield(opt,'tmean');       opt.tmean  = [];                    end;
if ~isfield(opt,'tval');        opt.tval = [];                      end;

dx = opt.dx;

time = (0:size(soma_pos,1)-1) * opt.tpf;
theta = zeros(size(time));
somadist = zeros(size(time));
bend = zeros(size(time));

%% compute Theta for Polar Coordinate & Bend parameter
st = track_time+1;
en = find(cellfun(@(x)~isempty(x),configs(:,1)),1,'last');
soma_pos_cycle = soma_pos(st:en,:);
soma_adj_cycle = soma_adj(st:en,:);
if any(all(soma_pos_cycle==0,2))
    error('No Soma position');
end
if any(all(soma_adj_cycle==0,2))
    error('No Soma position');
end
ds = sqrt(sum((soma_pos_cycle(:,1:2) - soma_adj_cycle(:,1:2)).^2,2));
[Dlarge, ~] = max(ds);
[Dsmall, small_pos] = min(ds);
d2theta = pi / (Dlarge - Dsmall);
dds = Dlarge - ds;
theta_cycle = dds * d2theta;
theta_cycle(small_pos:end) = 2*pi - theta_cycle(small_pos:end);
theta(track_time+1) = theta_cycle(1);
somadist(track_time+1) = ds(1)*opt.pixelsize;

bend(track_time+1) = sum(cellfun(@sum,configs(track_time+1,:)));

%% Init dendrite model
sizeI = [size(V,1), size(V,2)];
Zmax = max(length(configs{track_time+1,1}), length(configs{track_time+1,2}))*dx;
Zmin = -2;%length(configs{1,2})*dx;

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

%% create figure
f = figure('Visible','on','Position',[160,200,1050,785],'Color','w');
ha1 = axes('Units','normalized','Position',[0.6,0.05,0.35,0.4],'FontSize',16);
ha2 = axes('Units','normalized','Position',[0.6,0.55,0.35,0.4]);
ha3 = axes('Units','normalized','Position',[0.05,0.35,0.5,0.2]);
ha4 = axes('Units','normalized','Position',[0.05,0.08,0.5,0.2]);
ha5 = axes('Units','normalized','Position',[0.05,0.6,0.5,0.3]);

% plot polar coordinate
axes(ha1);
hradian = polar(ha1, theta_cycle(1),1,'.');
set(hradian, 'MarkerSize', 50);
set(ha1, 'Xdir', 'reverse');
hands = findall(f, 'parent', ha1, 'Type', 'text');
i = 10;
for hand = hands'
    if str2double(hand.String) > 0.1 && str2double(hand.String) < 1.1
        hand.String = [repmat(' ', 1, i-1) hand.String];
        i = i - 2;
    end
end
hradian_tlt = title(sprintf('Angle = %6.2f \\circ', theta_cycle(1)), 'FontSize', 16);
hold on;
htsem = patch(0,0,'r','EdgeAlpha', 0);
alpha(0.3)
htmean = polar(0,0);
set(htmean, 'color', 'k', 'linewidth',3)
hold off;
th = findall(gcf,'Type','text');
for i = 1:length(th),
    set(th(i),'FontSize',16)
end

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
axes(ha2);
hmodel_tlt = title(sprintf('Bending = %4.2f rad', bend(track_time+1)), 'FontSize', 16);
hmodel = patch(sf);
hmodel.FaceColor = 'red';
hmodel.EdgeColor = 'none';
daspect([1 1 1])
xlabel('X (\mum)', 'FontSize', 16);
ylabel('Y (\mum)', 'FontSize', 16);
zlabel('Z (\mum)', 'FontSize', 16);
planebox = [-Zmax Zmax yrange-somaY+[-10 10] Zmin Zmax] * opt.pixelsize;
axis([-Zmax Zmax yrange-somaY+[-10 10] Zmin Zmax] * opt.pixelsize)
zlim([Zmin Zmax] * opt.pixelsize)
view([15 36]);
set(ha2,'Ydir','reverse');
set(ha2,'Zdir','reverse');

hold on;
patch(planebox([1 2 2 1]), planebox([3 3 4 4]), planebox(5)*ones(1,4), zeros(1,4),'FaceAlpha',.1);
hold off;
box on
camlight 
lighting gouraud

% plot soma distance
axes(ha3);
htheta = plot(time, somadist, 'LineWidth', 3);
set(gca, 'FontSize', 16);
axis([([st en]-1)*opt.tpf 0 Dlarge*opt.pixelsize]);
hold on;
hdiv1 = plot([st-1 st-1]*opt.tpf, [0 Dlarge*opt.pixelsize], 'r', 'LineWidth', 2);
hold off;
ylabel('Soma Distance (\mum)', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);

% plot bend parameters
axes(ha4);
hbend = plot(time, bend, 'LineWidth', 3);
set(gca, 'FontSize', 16);
radmax = max(sum(cellfun(@sum, configs),2));
axis([([st en]-1)*opt.tpf 0 radmax]);
hold on;
hdiv2 = plot([st-1 st-1]*opt.tpf, [0 radmax], 'r', 'LineWidth', 2);
hold off;
ylabel('Bending (rad)', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);

% draw image and highlight dendrite
axes(ha5);
himg = imshow(repmat(V(:,:,track_time+1),1,1,3));
green = cat(3, zeros(sizeI), ones(sizeI), zeros(sizeI));
hold on;
hmask = imshow(green);
set(hmask, 'AlphaData', I2D_mask*.2);
ang = linspace(0,2*pi,100);
hsoma_pos = plot(soma_pos(track_time+1,1)+sin(ang)*soma_pos(track_time+1,3), ...
                 soma_pos(track_time+1,2)+cos(ang)*soma_pos(track_time+1,3), ...
                 'LineWidth', 2);
hsoma_adj = plot(soma_adj(track_time+1,1)+sin(ang)*soma_adj(track_time+1,3), ...
                 soma_adj(track_time+1,2)+cos(ang)*soma_adj(track_time+1,3), ...
                 'LineWidth', 2);hold off;
himg_tlt = title(sprintf('Time %.3f s, Frame %d', time(track_time+1), track_time), 'FontSize', 16);

% save to video
f.Visible = 'on';
if ~isempty(filename)
    writeVideo(outputVideo, getframe(f));
end

%% create animation
tval_ct = 1;
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
    set(himg_tlt, 'String', sprintf('Time %.3f s, Frame %d', time(t), t-1));
    set(hmask, 'AlphaData', I_mask_t*.2);
    set(hsoma_pos, 'XData', soma_pos(t,1)+sin(ang)*soma_pos(t,3));
    set(hsoma_pos, 'YData', soma_pos(t,2)+cos(ang)*soma_pos(t,3));
    set(hsoma_adj, 'XData', soma_adj(t,1)+sin(ang)*soma_adj(t,3));
    set(hsoma_adj, 'YData', soma_adj(t,2)+cos(ang)*soma_adj(t,3));
    
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

    % update theta
    somadist(t) = ds(t - track_time)*opt.pixelsize;
    theta(t) = theta_cycle(t - track_time);
    set(hradian_tlt, 'String', sprintf('Angle = %6.2f \\circ', theta(t)*180/pi));
    [xx, yy] = pol2cart(theta(t), 1);
    set(hradian, 'XData', xx);
    set(hradian, 'YData', yy);
    set(htheta, 'YData', somadist);
    set(hdiv1, 'XData', [t-1 t-1]*opt.tpf);
    
    % update bend
    bend(t) = sum(cellfun(@sum,configs(t,:)));
    set(hmodel_tlt, 'String', sprintf('Bending = %4.2f rad', bend(t)));
    set(hbend, 'YData', bend);
    set(hdiv2, 'XData', [t-1 t-1]*opt.tpf);
    
    newtval_ct = find(opt.tval < theta(t), 1, 'last');
    if isempty(newtval_ct), newtval_ct = 1; end
    tval_ct = max(tval_ct, newtval_ct);
    if t == size(V,3) || isempty(configs{t+1,1}), tval_ct = length(opt.tval); end
    
    if ~isempty(opt.tsem) && ~isempty(opt.tval) && ~isempty(opt.tmean)
        [x, y] = pol2cart([opt.tval(1:tval_ct) opt.tval(tval_ct:-1:1)], ...
                [opt.tmean(1:tval_ct) + opt.tsem(1:tval_ct) ...
                opt.tmean(tval_ct:-1:1) - opt.tsem(tval_ct:-1:1)]);
        set(htsem, 'XData', x);
        set(htsem, 'YData', y);
        [x, y] = pol2cart(opt.tval(1:tval_ct) ,opt.tmean(1:tval_ct));
        set(htmean, 'XData', x);
        set(htmean, 'YData', y);
    end

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

