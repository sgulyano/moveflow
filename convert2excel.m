function convert2excel( savedir, savename, gfp, num_img, initswctime, ddad_swc, ddae_swc, soma_cen, soma_rad, dd_plane, opt )
%CONVERT2EXCEL Summary of this function goes here
%   Detailed explanation goes here

if nargin < 10;  opt = struct();  end
if ~isfield(opt,'is_gui');              opt.is_gui      = true;         end
if ~isfield(opt,'ddad_left');           opt.ddad_left   = false;      	end
if ~isfield(opt,'deform_type');         opt.deform_type = 'prox-dist';  end
if ~isfield(opt,'directions');          opt.directions  = {};           end


savefile = fullfile(savedir, savename);

T = table((0:num_img)');
T.Properties.VariableNames{end} = 'ImageNumber';

bends = zeros(length(initswctime), 2, num_img+1);
masks = cell(length(initswctime), 2, num_img+1);
for i = 1:length(initswctime)
    opt.maxdep = 30;
    masks(i,:,:) = ffd_plot_both([], gfp, {ddad_swc{i}, ddae_swc{i}}, initswctime(i), ...
            {soma_cen{1,i}, soma_cen{2,i}}, dd_plane(i).configs, dd_plane(i).spline, opt);
    
    switch opt.deform_type
        case 'curvature'
            bend = cellfun(@(x)sum(abs((x))), results(num).configs);
        case 'prox-dist'
            bend = cellfun(@(x)sum(abs(cumsum(x))), dd_plane(i).configs);
    end
    
    if isempty(opt.directions)
        if opt.ddad_left
            bends(i,:,:) = bend';
        else
            bends(i,:,:) = bend(:,[2 1])';
        end
    else
        if opt.ddad_left
            bends(i,:,:) = ([opt.directions{1,i} -opt.directions{2,i}] .* bend)';
        else
            bends(i,:,:) = ([-opt.directions{1,i} opt.directions{2,i}] .* bend(:,[2 1]))';
        end
    end
    
    Td = table(soma_cen{1,i}(:,1), soma_cen{1,i}(:,2), squeeze(bends(i,1,:)));
    Td.Properties.VariableNames = {['XsomaddaD' num2str(i)], ['YsomaddaD' num2str(i)], ['CurvatureddaD' num2str(i)]};
    Te = table(soma_cen{2,i}(:,1), soma_cen{2,i}(:,2), squeeze(bends(i,2,:)));
    Te.Properties.VariableNames = {['XsomaddaE' num2str(i)], ['YsomaddaE' num2str(i)], ['CurvatureddaE' num2str(i)]};
    T = [T Td Te];
end
writetable(T, [savefile '.csv']);
disp(['Save track result to ' savefile '.csv ... Done']);

%% sanity check fitting neuron model
T = readtable([savefile '.csv']);
col1 = [235 90  235]/255; % ddaD color
col2 = [80  200 80 ]/255; % ddaE color
s = size(gfp);

somaX = (table2array(T(:,2:6:end)) + table2array(T(:,5:6:end))) / 2;
somaY = (table2array(T(:,3:6:end)) + table2array(T(:,6:6:end))) / 2;

outputVideo = VideoWriter([savefile '_' opt.deform_type '.mp4'], 'MPEG-4');
open(outputVideo);

f = figure(6); set(gcf, 'Pos', [100 100 1400 1000]);
subplot(3,1,1); himg = imshow(gfp(:,:,1), []); 
htitle = title('Frame 0', 'FontSize', 14);
hold on;
htxts = cell(1,length(initswctime));
for i = 1:length(initswctime)
    htxts{i} = text(somaX(1,i),somaY(1,i),num2str(i),'Color','red','FontSize',16);
    if somaX(1,i) == 0 && somaY(1,i) == 0
        htxts{i}.Visible = 'off';
    else
        htxts{i}.Visible = 'on';
    end
end

hmask1 = imshow(repmat(reshape(col1, [1 1 3]), s([1 2])));
I2D_mask1 = zeros(s([1 2]));
for i = 1:length(initswctime)
    if ~isempty(masks{i,1,1})
        I2D_mask1 = I2D_mask1 | masks{i,1,1};
    end
end
set(hmask1, 'AlphaData', I2D_mask1*.2);
hmask2 = imshow(repmat(reshape(col2, [1 1 3]), s([1 2])));
I2D_mask2 = zeros(s([1 2]));
for i = 1:length(initswctime)
    if ~isempty(masks{i,2,1})
        I2D_mask2 = I2D_mask2 | masks{i,2,1};
    end
end
set(hmask2, 'AlphaData', I2D_mask2*.2);

maxval = max(bends(:));
minval = min(bends(:));
subplot(3,1,2); plot(0:num_img, table2array(T(:,4:6:end)));
hold on;
hcur1 = plot([0 0], [minval maxval], 'r');
hold off;
axis tight
switch opt.deform_type
    case 'curvature'
        title('Sum Abs Curvature of ddaD', 'FontSize', 14);
    case 'prox-dist'
        title('Acc Curvature Proximal->Distal ddaD', 'FontSize', 14);
end

subplot(3,1,3); plot(0:num_img, table2array(T(:,7:6:end)));
hold on;
hcur2 = plot([0 0], [minval maxval], 'r');
hold off;
axis tight
switch opt.deform_type
    case 'curvature'
        title('Sum Abs Curvature of ddaE', 'FontSize', 14);
    case 'prox-dist'
        title('Acc Curvature Proximal->Distal ddaE', 'FontSize', 14);
end

[lgd, hobj] = legend(cellfun(@num2str, num2cell(1:length(initswctime))'));
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',5);
ht = findobj(hobj,'type','text');
set(ht,'FontSize',14);
lgd.FontSize = 10;

writeVideo(outputVideo, getframe(f));

for t = 2:num_img+1
    set(himg, 'CData', gfp(:,:,t));
    set(htitle, 'String', ['Frame ' num2str(t-1)]);
    set(hcur1, 'XData', [t t]-1);
    set(hcur2, 'XData', [t t]-1);
    
    for i = 1:length(initswctime)
        set(htxts{i}, 'Position', [somaX(t,i),somaY(t,i),0]);
        if somaX(t,i) == 0 && somaY(t,i) == 0
            htxts{i}.Visible = 'off';
        else
            htxts{i}.Visible = 'on';
        end
    end
    
    I2D_mask1 = zeros(s([1 2]));
    for i = 1:length(initswctime)
        if ~isempty(masks{i,1,t})
            I2D_mask1 = I2D_mask1 | masks{i,1,t};
        end
    end
    set(hmask1, 'AlphaData', I2D_mask1*.2);
    
    I2D_mask2 = zeros(s([1 2]));
    for i = 1:length(initswctime)
        if ~isempty(masks{i,2,t})
            I2D_mask2 = I2D_mask2 | masks{i,2,t};
        end
    end
    set(hmask2, 'AlphaData', I2D_mask2*.2);
    
    drawnow;
    writeVideo(outputVideo, getframe(f));
end

disp(['Save video result to ' savefile '_' opt.deform_type '.mp4 ... Done']);
if opt.is_gui
    uiwait(msgbox('Done'));
else
    pause;
end
delete(f);
end

