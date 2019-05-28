close all; clear all; clc;

addpath(genpath('toolbox'));
addpath(genpath('UGM'))

functionname = 'trackGFPRFPXYZTdata.m';
functiondir = which(functionname);
functiondir = functiondir(1:end-length(functionname));

%% user parameters
DEBUG = true;
track_opt.DEBUG = DEBUG;
SAVE_VIDEO = true;
dataset = 8;
FULLFLOW = false;
switch dataset
    case 8
        gfpfile = '~/Desktop/LipingsData/20171228L2/20171228L2_z%02d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 635;                                  % number of frames (start from 0 to num_img)
        savefile = '20171228L2/20171228L2';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = [1 397 462 511];                   % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = false;                              % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 9
        gfpfile = '~/Desktop/LipingsData/20171228L1_B/20171228L1_B_z%02d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 12;                                 % number of image slices per stack
        end_slice = 12;                                 % image slice to stop reading
        num_img = 201;                                  % number of frames (start from 0 to num_img)
        savefile = '20171228L1_B/20171228L1_B';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 1;                                 % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = false;                              % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 10
        gfpfile = '~/Desktop/LipingsData/20171228L5/20171228L5_z%02d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 328;                                  % number of frames (start from 0 to num_img)
        savefile = '20171228L5/20171228L5';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = [1 45 132 255];                    % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = true;                               % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 11
        gfpfile = '~/Desktop/LipingsData/20171229L2/20171229L2_z%02d_t%02d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 92;                                   % number of frames (start from 0 to num_img)
        savefile = '20171229L2/20171229L2';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 1;                                 % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = true;                               % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 12
        gfpfile = '~/Desktop/LipingsData/20171229L8-B/20171229L8_B_good_z%02d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 205;                                  % number of frames (start from 0 to num_img)
        savefile = '20171229L8_B/20171229L8_B';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 60;                                % frame picked for manual tracing
        pickslice = 5;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = true;                               % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 13
        gfpfile = '~/Desktop/LipingsData/20171229L12_B/20171229L12_B_good_z%02d_t%02d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 71;                                  % number of frames (start from 0 to num_img)
        savefile = '20171229L12_B/20171229L12_B';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 1;                                 % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = true;                               % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 14
        gfpfile = '~/Desktop/LipingsData/20171229L13/20171229L13_F_z%02d_t%02d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 11;                                 % number of image slices per stack
        end_slice = 11;                                 % image slice to stop reading
        num_img = 106;                                  % number of frames (start from 0 to num_img)
        savefile = '20171229L13/20171229L13';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = [1 35];                            % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = true;                               % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
end

%% read img data
disp('Reading GFP and RFP images');
gfpinfo = imfinfo(sprintf(gfpfile,0,0));
gfp = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
gfp1sl = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
if ~isempty(rfpfile)
    rfpinfo = imfinfo(sprintf(rfpfile,0,0));
    if gfpinfo.Height ~= rfpinfo.Height || gfpinfo.Width ~= rfpinfo.Width,
        error('GFP and RFP size mismatch');
    end
    rfp = zeros(rfpinfo.Height, rfpinfo.Width, num_img+1, 'uint8');
end

for t = 0:num_img
    V = zeros(gfpinfo.Height, gfpinfo.Width, num_slice+1, 'uint8');
    for z = 0:num_slice
        if switchTZ
            X = imread( sprintf(gfpfile,z,t) );
        else
            X = imread( sprintf(gfpfile,t,z) );
        end
        if normalized
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        end
        if z == pickslice
            gfp1sl(:,:,t+1) = X;
        end
        if adj_slicewise
            n = 3;
            Idouble = im2double(X);
            avg = mean2(Idouble);
            sigma = std2(Idouble);
            X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
        end
        V(:,:,z+1) = X;
        
        if any(t == pickframes) && DEBUG
            pickframedir = [savefile '_ref_T' num2str(t)];
            if ~exist(pickframedir, 'dir') && DEBUG
                mkdir(pickframedir);
            end
            imwrite(X, sprintf([pickframedir '/frame%02d.tif'], z));
        end
    end
    gfp(:,:,t+1) = max(V(:,:,start_slice+1:end_slice+1),[],3);
    
    if ~isempty(rfpfile)
        V = zeros(rfpinfo.Height, rfpinfo.Width, num_slice+1);
        for z = 0:num_slice
            V(:,:,z+1) = imread( sprintf(rfpfile,t,z) );
        end
        rfp(:,:,t+1) = max(V,[],3);
    end
    fprintf('.');
    if mod(t+1,100) == 0, fprintf('\n'); end;
end
fprintf('\n');

if DEBUG
    figure(1); clf; imshow3D(gfp); title(savefile);
%     imgflow2vdo( [savefile '_gfp.mp4'], gfp );
    load([savefile '_align_gfp.mat'], 'gfp');
    figure(2); clf; imshow3D(gfp); title('align');
%     imgflow2vdo( [savefile '_align_gfp.mp4'], gfp );
%     keyboard;
    delete(1);
    delete(2);
end

%% get swc list
swclist = dir([savefile '*.swc']);
initswctime = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
initswctime = cellfun(@(x)str2double(x{1}), initswctime);

if exist([savefile '_somaboth.mat'], 'file')
    load([savefile '_somaboth.mat']);
else
    %% track soma using hough transform
    % track soma for each trace
    disp('Tracking soma... started');
    centers = cell(num_img+1,1);
    radii = cell(num_img+1,1);
    for t = 1:num_img+1
        [centers{t}, radii{t}] = imfindcircles(gfp1sl(:,:,t), soma_opt.radii, ...
                'ObjectPolarity', 'bright', ...
                'Sensitivity', soma_opt.sensitiv, ...
                'EdgeThreshold',soma_opt.edge_thr, ...
                'Method','TwoStage');
    end
    %%
    soma_pos = cell(length(initswctime),1);
    for num = 1:length(initswctime)
        track_time = initswctime(num);
        swcs = read_swc_file( fullfile(fileparts(savefile), swclist(num).name) );
        swcs(:,3:5) = swcs(:,3:5)+1;
        
        rootidx = [find(swcs(:,7)==-1); size(swcs,1)+1];
        
        tokens = strsplit(swclist(num).name, '_');
        
        neuron = arrayfun(@(K)struct('type', tokens(end-1), 'init_time', track_time, ...
                'cen', zeros(num_img+1, 2), 'rad', zeros(num_img+1, 1), 'swc', []), ...
                1:length(rootidx)-1, 'UniformOutput',true);
        
        for i = 1:length(rootidx)-1
            swc = swcs(rootidx(i):rootidx(i+1)-1,:);
            swc(:,[1 7]) = swc(:,[1 7]) - rootidx(i) + 1;
            swc(1,7) = -1;
            neuron(i).swc = swc;
            
            [~, pos] = min(sum(bsxfun(@minus, centers{track_time+1}, swc(1,3:4)).^2,2));
            neuron(i).cen(track_time+1,:) = centers{track_time+1}(pos,:);
            neuron(i).rad(track_time+1) = radii{track_time+1}(pos);

            figure(5);
            EV_plot_img(gfp1sl(:,:,track_time+1), swc); 
            hold on;
            EV_plot_img([], swc);
            title(num2str(track_time));
            viscircles(centers{track_time+1}, radii{track_time+1}, 'EdgeColor', 'b');
            viscircles(neuron(i).cen(track_time+1,:), neuron(i).rad(track_time+1), 'EdgeColor', 'r');
            drawnow;
            title([neuron(i).type '-t' num2str(track_time)]);
            pause;
            % track forward
            [neuron(i).cen, neuron(i).rad] = soma_manual( gfp1sl, neuron(i).cen, neuron(i).rad, track_time+1, 1, centers, radii );
        end
        soma_pos{num} = neuron;
    end
    soma_pos = [soma_pos{:}];
    save([savefile '_somaboth.mat'], 'soma_pos');
    disp(['Tracking soma... completed : ' savefile '_somaboth.mat']);
end

%% sanity check soma
col = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
for t = 1:num_img+1
    figure(5); imshow(gfp(:,:,t), []); title(t-1);
    for ii = 1:length(soma_pos)
        if ~all(soma_pos(ii).cen(t,:) == 0)
            viscircles(soma_pos(ii).cen(t,:), soma_pos(ii).rad(t), 'EdgeColor', col(mod(ii-1,7)+1,:));
        end
    end
    drawnow;
end

%% use GUI to fit the folding neurons
if exist([savefile '_model_both.mat'], 'file')
    load([savefile '_model_both.mat'], 'results');
else
    %%
    results = struct([]);
    ddad_idx = cellfun(@(x)strcmp(x,'ddaD'), {soma_pos.type});
    soma_ddad = soma_pos(ddad_idx);
    soma_ddae = soma_pos(~ddad_idx);
    for num = 1:length(soma_ddad)
        track_time = soma_ddad(num).init_time;
        swc_ddad = soma_ddad(num).swc;

        soma_ddae_track_time = soma_ddae([soma_ddae.init_time] == track_time);
        ddae_pos = arrayfun(@(x)x.swc(1,3:5), soma_ddae_track_time, 'UniformOutput', false);
        ds = bsxfun(@minus, vertcat(ddae_pos{:}), swc_ddad(1,3:5));
        [~, pos] = min(sum(ds.^2,2));

        swc_ddae = soma_ddae_track_time(pos).swc;
        
        % fit neuron
        [configs,leftpos,rightpos,spline] = ffd_gui_both(gfp, {swc_ddad, swc_ddae}, ...
                track_time, {soma_ddad(num).cen, soma_ddae_track_time(pos).cen});
        videoname = sprintf([savefile '_model_both%d.mp4'], num);
        ffd_plot_both(videoname, gfp, {swc_ddad, swc_ddae}, track_time, {soma_ddad(num).cen, soma_ddae_track_time(pos).cen}, configs, spline);
        results(num).configs = configs;
        results(num).leftpos = leftpos;
        results(num).rightpos = rightpos;
        results(num).spline = spline;
    end
    save([savefile '_model_both.mat'], 'results');
    disp(['Save tracking to ' savefile '_model_both.mat Done'])
end

%% convert to Excel format
% load([savefile '_model_both.mat'], 'results');
load ddaeforward_back.mat
load thetain.mat
opt.tsem = ddaeforward_back_sem;
opt.tmean = ddaeforward_back_mean;
opt.tval = thetain;

T = table((0:num_img)');
T.Properties.VariableNames{end} = 'ImageNumber';

ddad_idx = cellfun(@(x)strcmp(x,'ddaD'), {soma_pos.type});
soma_ddad = soma_pos(ddad_idx);
soma_ddae = soma_pos(~ddad_idx);
bends = zeros(length(initswctime), 2, num_img+1);
masks = cell(length(initswctime), 2, num_img+1);
for num = 1:length(initswctime)
    opt.maxdep = 30;
    track_time = soma_ddad(num).init_time;
    swc_ddad = soma_ddad(num).swc;

    soma_ddae_track_time = soma_ddae([soma_ddae.init_time] == track_time);
    ddae_pos = arrayfun(@(x)x.swc(1,3:5), soma_ddae_track_time, 'UniformOutput', false);
    ds = bsxfun(@minus, vertcat(ddae_pos{:}), swc_ddad(1,3:5));
    [~, pos] =min(sum(ds.^2,2));

    swc_ddae = soma_ddae_track_time(pos).swc;

%     videoname = sprintf([savefile '_model_both%d.mp4'], num);
    masks(num,:,:) = ffd_plot_both([], gfp, {swc_ddad, swc_ddae}, track_time, {soma_ddad(num).cen, soma_ddae_track_time(pos).cen}, results(num).configs, results(num).spline, opt);
    
%     videoname = sprintf([savefile '_model_both_stat%d.mp4'], num);
%     theta_ffd_plot_both(videoname, gfp, {swc_ddad, swc_ddae}, track_time, {soma_ddad(num).cen, soma_ddae(pos).cen}, soma_adj{num}, results(num).configs, results(num).spline, opt);


%     bend = cellfun(@(x)sum(abs((x))), results(num).configs);
    bend = cellfun(@(x)sum(abs(cumsum(x))), results(num).configs);
    
    if ddaD_Left
        bends(num,:,:) = bend';
    else
        bends(num,:,:) = bend(:,[2 1])';
    end
    
    Td = table(soma_ddad(num).cen(:,1), soma_ddad(num).cen(:,2), squeeze(bends(num,1,:)));
    Td.Properties.VariableNames = {['XsomaddaD' num2str(num)], ['YsomaddaD' num2str(num)], ['CurvatureddaD' num2str(num)]};
    Te = table(soma_ddae_track_time(pos).cen(:,1), soma_ddae_track_time(pos).cen(:,2), squeeze(bends(num,2,:)));
    Te.Properties.VariableNames = {['XsomaddaE' num2str(num)], ['YsomaddaE' num2str(num)], ['CurvatureddaE' num2str(num)]};
    T = [T Td Te];
end
writetable(T, [savefile '_acmcurv.csv']);
disp(['Save track result to ' savefile '_acmcurv.csv ... Done']);

%% sanity check fitting neuron model
T = readtable([savefile '_acmcurv.csv']);
col1 = [235 90  235]/255; % ddaD color
col2 = [80  200 80 ]/255; % ddaE color
s = size(gfp);

somaX = (table2array(T(:,2:6:end)) + table2array(T(:,5:6:end))) / 2;
somaY = (table2array(T(:,3:6:end)) + table2array(T(:,6:6:end))) / 2;

outputVideo = VideoWriter([savefile '_acmcurv.mp4'], 'MPEG-4');
open(outputVideo);

f = figure(6); set(gcf, 'Pos', [100 100 1400 1000]);
subplot(3,1,1); himg = imshow(gfp(:,:,1), []); 
htitle = title('Frame 0', 'FontSize', 14);
hold on;
htxts = cell(1,length(initswctime));
for i = 1:length(initswctime)
    htxts{i} = text(somaX(1,i),somaY(1,i),num2str(i),'Color','red','FontSize',16);
    if somaX(1,i) == 0 && somaY(1,i) == 0,
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
subplot(3,1,2); plot(0:num_img, table2array(T(:,4:6:end)));
hold on;
hcur1 = plot([0 0], [0 maxval], 'r');
hold off;
axis tight
title('ddaD', 'FontSize', 14);

subplot(3,1,3); plot(0:num_img, table2array(T(:,7:6:end)));
hold on;
hcur2 = plot([0 0], [0 maxval], 'r');
hold off;
axis tight
title('ddaE', 'FontSize', 14);

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
        if somaX(t,i) == 0 && somaY(t,i) == 0,
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
close(outputVideo);
disp(['Save video result to ' savefile '_acmcurv.mp4 ... Done']);
pause;
delete(f);