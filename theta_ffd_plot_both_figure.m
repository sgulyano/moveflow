function theta_ffd_plot_both( filename, V, swcs, track_time, soma_pos, soma_adj, configs, spline, opt )
%THETA_FFD_PLOT_BOTH draw results tracking in 3D over Depth using FFD and 
%save to video

if nargin < 9;  opt = struct();  end;
if ~isfield(opt,'dx');          opt.dx = 4;                         end;
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end;
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end;
if ~isfield(opt,'tsem');        opt.tsem   = [];                    end;
if ~isfield(opt,'tmean');       opt.tmean  = [];                    end;
if ~isfield(opt,'tval');        opt.tval = [];                      end
if ~isfield(opt,'maxdep');      opt.maxdep    = 27;                 end;
if ~isfield(opt,'zflip');       opt.zflip     = [];                 end;

col1 = [235 90  235]/255; % ddaD color
col2 = [80  200 80 ]/255; % ddaE color

dx = opt.dx;

time = (0:size(soma_pos{1},1)-1) * opt.tpf;
theta = zeros(size(time));
somadist = zeros(size(time));

%% compute Theta for Polar Coordinate & Bend parameter
st = track_time+1;
en = find(cellfun(@(x)~isempty(x),configs(:,1)),1,'last');
ddad_pos_cycle = soma_pos{1}(st:en,:);
ddae_pos_cycle = soma_pos{2}(st:en,:);
soma_adj_cycle = soma_adj(st:en,:);
if any(all(ddad_pos_cycle==0,2)) || any(all(ddae_pos_cycle==0,2))
    error('No ddaD/ddaE soma position');
end
if any(all(soma_adj_cycle==0,2))
    error('No adjacent soma position');
end
ds = sqrt(sum((ddad_pos_cycle(:,1:2) - soma_adj_cycle(:,1:2)).^2,2));
[Dlarge, ~] = max(ds);
[Dsmall, small_pos] = min(ds);
d2theta = pi / (Dlarge - Dsmall);
dds = Dlarge - ds;
theta_cycle = dds * d2theta;
theta_cycle(small_pos:end) = 2*pi - theta_cycle(small_pos:end);
theta(track_time+1) = theta_cycle(1);
somadist(track_time+1:track_time+length(ds)) = ds*opt.pixelsize;
bend = cellfun(@(x)sum(abs(x)),configs);

%% Init dendrite model
sizeI = [size(V,1), size(V,2)];
Xmax = max(length(configs{track_time+1,1}), length(configs{track_time+1,2}))*dx;
Zmax = opt.maxdep*dx;

% create baseline mask
I_mask1 = swc2pixel( swcs{1}, sizeI );
I_mask2 = swc2pixel( swcs{2}, sizeI );
I2D_mask = imdilate(I_mask1 | I_mask2, ones(3));
I2D_mask1 = imdilate(I_mask1, ones(3));
I2D_mask2 = imdilate(I_mask2, ones(3));
% get bounding box
[ii, jj] = find(I2D_mask);
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
I_mask1_crop = I_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
I_mask2_crop = I_mask2(yrange(1):yrange(2),xrange(1):xrange(2));
I2D_mask1_crop = I2D_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
I2D_mask2_crop = I2D_mask2(yrange(1):yrange(2),xrange(1):xrange(2));
% locate soma
ddad_somaX = round(swcs{1}(1,3));
ddad_somaY = round(swcs{1}(1,4));
ddae_somaX = round(swcs{2}(1,3));
ddae_somaY = round(swcs{2}(1,4));
somaX = round((ddad_somaX + ddae_somaX)/2);
somaY = round((ddad_somaY + ddae_somaY)/2);

% init FFD
if isempty(opt.zflip)
    [Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx);
else
    [Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx, opt.zflip{track_time+1});
end

% open video
if ~isempty(filename)
    outputVideo = VideoWriter(filename, 'MPEG-4');
    open(outputVideo);
end

%% create figure
f = figure('Visible','on','Position',[160,200,1050,785],'Color','w');
ha1 = axes('Units','normalized','Position',[0.6,0.05,0.35,0.4],'FontSize',16);
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
ddad_somaZ = cz(ddad_somaY - yrange(1) + 1, ddad_somaX - xrange(1) + 1);
ddae_somaZ = cz(ddae_somaY - yrange(1) + 1, ddae_somaX - xrange(1) + 1);
cx = cx(:); cy = cy(:); cz = cz(:);

% create 3D surface of dendrite morphology
sizeI3D = [sizeI 1+range(cz)];
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(I_mask1_crop), cx(I_mask1_crop), cz(I_mask1_crop))) = 1;
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
[xx3, yy3, zz3] = meshgrid(1:size(I3D,2), 1:size(I3D,1), 1:size(I3D,3));
I3D((xx3 - ddad_somaX - 2).^2 + (yy3 - ddad_somaY - 2).^2 + (zz3 - ddad_somaZ - 7).^2 < swcs{1}(1,6)^2) = 1;
sf1 = isosurface(I3D, .5);
% relocate soma to (0,0,1)
% sf1.vertices = bsxfun(@minus, sf1.vertices, [somaX, somaY, 0]) * opt.pixelsize;

% create 3D surface of ddaE
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(I_mask2_crop), cx(I_mask2_crop), cz(I_mask2_crop))) = 1;
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
I3D((xx3 - ddae_somaX - 2).^2 + (yy3 - ddae_somaY - 2).^2 + (zz3 - ddae_somaZ - 7).^2 < swcs{2}(1,6)^2) = 1;
sf2 = isosurface(I3D, .5);
% relocate soma to (0,0,1)
% sf2.vertices = bsxfun(@minus, sf2.vertices, [somaX, somaY, 0]) * opt.pixelsize;

I3D = zeros(sizeI3D);
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D((xx3 - soma_adj(track_time+1,1) - 2).^2 + (yy3 - soma_adj(track_time+1,2) - 2).^2 + (zz3 - soma_adj(track_time+1,3) - 7).^2 < swcs{2}(1,6)^2) = 1;
sf3 = isosurface(I3D, .5);


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
plot(time, bend(:,2), 'LineWidth', 3, 'Color', col1);
hold on;
plot(time, bend(:,1), 'LineWidth', 3, 'Color', col2);
set(gca, 'FontSize', 16);
radmax = max(max(cellfun(@(x)sum(abs(x)), configs)));
axis([([st en]-1)*opt.tpf 0 radmax]);
hdiv2 = plot([st-1 st-1]*opt.tpf, [0 radmax], 'r', 'LineWidth', 2);
hold off;
ylabel('Bending (rad)', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
legend({'ddaD', 'ddaE'}, 'Location', 'northwest');

% draw image and highlight dendrite
axes(ha5);
[xx, yy, zz] = meshgrid(1:size(V,2), 1:size(V,1), Zmax * opt.pixelsize);
himg = surf(xx, yy, zz, V(:,:,track_time+1),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct',...
        'SpecularStrength',0,...
        'AmbientStrength',1);
colormap(gray(256));
axis tight
hold on;
hmodel1 = patch(sf1);
hmodel1.FaceColor = col1;
hmodel1.EdgeColor = 'none';
hmodel1.AmbientStrength = 0.5;

hmodel2 = patch(sf2);
hmodel2.FaceColor = col2;
hmodel2.EdgeColor = 'none';
hmodel2.AmbientStrength = 0.5;

hmodel3 = patch(sf3);
hmodel3.FaceColor = 'Cyan';
hmodel3.EdgeColor = 'none';
hmodel3.AmbientStrength = 0.5;
hold off;

daspect([1 1 1])
xlabel('X (\mum)', 'FontSize', 16);
ylabel('Y (\mum)', 'FontSize', 16);
zlabel('Z (\mum)', 'FontSize', 16);
zlim([0 Zmax] * opt.pixelsize)
view([0 90]);
box off
axis off
camproj('perspective')
hl = camlight('headlight');
lighting gouraud

set(ha5,'Ydir','reverse');
set(ha5,'Zdir','reverse');
camlight(hl,'left')

camorbit(20,0,'data',[1 0 0]);
camdolly(0,0,0.7);
% keyboard;
% 
himg_tlt = title(sprintf('Time %.3f s, Frame %d', time(track_time+1), track_time), 'FontSize', 16);
% axes(ha6)
% plot(theta_cycle);

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
        
    % update image and dendrite highlight
    set(himg, 'CData', V(:,:,t));
    set(himg_tlt, 'String', sprintf('Time %.3f s, Frame %d', time(t), t-1));
    
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
    % relocate soma to (0,0,1)
    hmodel1.Vertices = sf1.vertices;%bsxfun(@minus, sf1.vertices, [newsomaX, newsomaY, 0]) * opt.pixelsize;
    hmodel1.Faces = sf1.faces;
    % update 3D dendrite ddaE
    I3D = zeros(sizeI3D);
    I3D(sub2ind(sizeI3D, cy(I_mask2_crop(:)&idxIn), cx(I_mask2_crop(:)&idxIn), cz(I_mask2_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    I3D((xx3 - ddaeX - 2).^2 + (yy3 - ddaeY - 2).^2 + (zz3 - ddaeZ - 7).^2 < swcs{2}(1,6)^2) = 1;
    sf2 = isosurface(I3D, .5);
    % relocate soma to (0,0,1)
    hmodel2.Vertices = sf2.vertices;%bsxfun(@minus, sf2.vertices, [newsomaX, newsomaY, 0]) * opt.pixelsize;
    hmodel2.Faces = sf2.faces;
    
    I3D = zeros(sizeI3D);
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D((xx3 - soma_adj(t,1) - 2).^2 + (yy3 - soma_adj(t,2) - 2).^2 + (zz3 - soma_adj(t,3) - 7).^2 < swcs{2}(1,6)^2) = 1;
    sf3 = isosurface(I3D, .5);
    hmodel3.Vertices = sf3.vertices;
    hmodel3.Faces = sf3.faces;
    

    % update theta
    theta(t) = theta_cycle(t - track_time);
    set(hradian_tlt, 'String', sprintf('Angle = %6.2f \\circ', theta(t)*180/pi));
    [xx, yy] = pol2cart(theta(t), 1);
    set(hradian, 'XData', xx);
    set(hradian, 'YData', yy);
    set(hdiv1, 'XData', [t-1 t-1]*opt.tpf);
    
    % update bend
%     set(hmodel_tlt, 'String', sprintf('Bend: ddaE = %4.2f rad, ddaD = %4.2f rad', bend(t,:)));
    set(hdiv2, 'XData', [t-1 t-1]*opt.tpf);
    
    newtval_ct = find(opt.tval < theta(t), 1, 'last');
    if isempty(newtval_ct), newtval_ct = 1; end
    tval_ct = max(tval_ct, newtval_ct);
    if t == size(V,3) || isempty(configs{t+1,1}), tval_ct = length(opt.tval); end
    
    [x, y] = pol2cart([opt.tval(1:tval_ct) opt.tval(tval_ct:-1:1)], ...
            [opt.tmean(1:tval_ct) + opt.tsem(1:tval_ct) ...
            opt.tmean(tval_ct:-1:1) - opt.tsem(tval_ct:-1:1)]);
    set(htsem, 'XData', x);
    set(htsem, 'YData', y);
    [x, y] = pol2cart(opt.tval(1:tval_ct) ,opt.tmean(1:tval_ct));
    set(htmean, 'XData', x);
    set(htmean, 'YData', y);

    % save to video
    if t > 290
        pause;%(0.1);
    else
        pause(0.1);
    end
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

