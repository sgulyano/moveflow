% close all; clear all; clc;

addpath(genpath('Full_Flow_Source_Code'))
addpath(genpath('toolbox'));
addpath(genpath('UGM'))

functionname = 'trackGFPRFPXYZTdata.m';
functiondir = which(functionname);
functiondir = functiondir(1:end-length(functionname));

%% user parameters
DEBUG = true;
track_opt.DEBUG = DEBUG;
SAVE_VIDEO = true;
dataset = 7;
FULLFLOW = false;
switch dataset
    case 1
        gfpfile = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/larva3_2_z%d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        use_normxcorr = false;                          % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 8;                                  % number of image slices per stack
        end_slice = num_slice;                          % image slice to stop reading
        num_img = 409;                                  % number of frames (start from 0 to num_img)
        savefile = 'larva3_frames/larva3';
        gfpadjrange = [0.02,0.25];                      % intensity adjusting threshold for Ca images
        pickframes = [198, 201, 229, 252, 340, 342];    % frame picked for manual tracing
        pickslice = 4;                                  % image slice for tracking soma
        use_medfilt = false;                            % use median filter in soma detection
        soma_opt = struct();                            % options for soma detector
    case 2
        gfpfile = '~/Desktop/LipingsData/1121larva1_1/1122larva1_1good_z%d_t%03d.tif';
        rfpfile = '';
        switchTZ = true;
        start_slice = 0;
        num_slice = 8;
        end_slice = num_slice;                          % image slice to stop reading
        num_img = 751;
        savefile = '1122larva1/1122larva1';
        gfpadjrange = [0.0,0.1];
        pickframes = [43, 47, 184, 285, 445, 555, 660];
        pickslice = 4;                                  % image slice for tracking soma
        soma_opt = struct();                            % options for soma detector
    case 3
        gfpfile = '~/Desktop/LipingsData/1202DLa2_1good_subset/1202DLa2_1good_Subset_z%d_t%03d.tif';
        rfpfile = '';
        switchTZ = true;
        normalized = false;                             % normalized image slice
        start_slice = 3;
        num_slice = 7;
        end_slice = num_slice;                          % image slice to stop reading
        num_img = 899;
        savefile = '1202DLa2/1202DLa2';
        gfpadjrange = [0.17,0.7];
        pickslice = 4;                                  % image slice for tracking soma
        soma_opt = struct();                            % options for soma detector
    case 4
        gfpfile = '~/Desktop/LipingsData/GFPRFPXYZTdata/larva4S14Zgreen/Larva 4_Series014_Crop001_t%03d_z%d_ch00.tif';
        rfpfile = '~/Desktop/LipingsData/GFPRFPXYZTdata/Larva4S14redZ/Larva 4_Series014_Crop002_t%03d_z%d_ch00.tif';
        switchTZ = false;
        normalized = false;                             % normalized image slice
        use_normxcorr = false;                          % use Normalized Cross-Correlation to track soma
        start_slice = 0;
        num_slice = 2;      % number of image slices per stack
        end_slice = num_slice;                          % image slice to stop reading
        num_img = 292;      % number of frames (start from 0 to num_img)
        savefile = 'Larva4s014/Larva4s014';
        gfpadjrange = [0.04 0.5]; % for FullFlow use [0.04 0.3]
        rfpadjrange = [0 .5];
        pickframes = 1;     % frame picked for manual tracing
        pickslice = 1;      % image slice for tracking soma
        use_medfilt = '3D';
        soma_opt.radii = [25 50];
        soma_opt.sensitiv = .8;
    case 5
        gfpfile = '~/Desktop/LipingsData/1122larva5_3_Subset for Sarun/1122larva5_3_Subset_z%d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 6;                                  % number of image slices per stack
        end_slice = 5;                                  % image slice to stop reading
        num_img = 256;                                  % number of frames (start from 0 to num_img)
        savefile = '1122larva5_3/1122larva5_3';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 152;                               % frame picked for manual tracing
        pickslice = 3;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        soma_opt = struct();                            % options for soma detector
        soma_opt.radii = [15 30];
        soma_opt.sensitiv = 0.6;
        soma_opt.edge_thr = 0.3;
    case 6
        gfpfile = '~/Desktop/LipingsData/20170616L3/20170616L3_z%02d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 10;                                 % number of image slices per stack
        end_slice = 10;                                 % image slice to stop reading
        num_img = 711;                                  % number of frames (start from 0 to num_img)
        savefile = '20170616L3/20170616L3';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 72;                                % frame picked for manual tracing
        motion_end_frames = 111;                        % frame that motion end correspond to pickframes
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        soma_opt = struct();                            % options for soma detector
        soma_opt.radii = [15 30];
        soma_opt.sensitiv = 0.6;
        soma_opt.edge_thr = 0.3;
        adj_soma = [29, 86];
    case 7
        gfpfile = 'C:\Users\sguly\Documents\MATLAB\DanData\20170616L7_align\20170616L7_align_z%d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 9;                                  % number of image slices per stack
        end_slice = 9;                                  % image slice to stop reading
        num_img = 791;                                  % number of frames (start from 0 to num_img)
        savefile = '20170616L7/20170616L7';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 48;                                % frame picked for manual tracing
        motion_end_frames = 150;                        % frame that motion end correspond to pickframes
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
        adj_soma = [140, 71, 5];
end

%% read img data
disp('Reading GFP and RFP images');
[di, ff, ext] = fileparts(gfpfile);
gfpinfo = imfinfo(fullfile(di, [sprintf(ff,0,0) ext]));
gfp = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
gfp1sl = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
if ~isempty(rfpfile)
    [rdi, rff, rext] = fileparts(gfpfile);
    rfpinfo = imfinfo(fullfile(rdi, [sprintf(rff,0,0) rext]));
    if gfpinfo.Height ~= rfpinfo.Height || gfpinfo.Width ~= rfpinfo.Width
        error('GFP and RFP size mismatch');
    end
    rfp = zeros(rfpinfo.Height, rfpinfo.Width, num_img+1, 'uint8');
end
for t = 0:num_img
    V = zeros(gfpinfo.Height, gfpinfo.Width, num_slice+1, 'uint8');
    for z = 0:num_slice
        if switchTZ
            X = imread( fullfile(di, [sprintf(ff,z,t) ext]) );
        else
            X = imread( fullfile(di, [sprintf(ff,t,z) ext]) );
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

%% get swc list
swclist = dir([savefile '_ddad_t*.swc']);
initswctime = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
initswctime = cellfun(@(x)str2double(x{1}), initswctime);

if exist([savefile '_somaboth.mat'], 'file')
    load([savefile '_somaboth.mat']);
else
    %% track soma using hough transform
    % track soma for each trace
    disp('Tracking soma... started');
    soma_pos = cell(length(initswctime),2);
    soma_adj = cell(length(initswctime),1);
    for num = 1%:length(initswctime)
        track_time = initswctime(num);
        swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
        swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
        swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
        swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;
        
        [centers, radii, metric] = imfindcircles(gfp1sl(:,:,track_time+1), soma_opt.radii, ...
                'ObjectPolarity', 'bright', ...
                'Sensitivity', soma_opt.sensitiv, ...
                'EdgeThreshold',soma_opt.edge_thr, ...
                'Method','TwoStage');
            
        [~, pos] = min(sum(bsxfun(@minus, centers, swc_ddad(1,3:4)).^2,2));
        ddad_cen = centers(pos,:);
        ddad_rad = radii(pos);
        
        [~, pos] = min(sum(bsxfun(@minus, centers, swc_ddae(1,3:4)).^2,2));
        ddae_cen = centers(pos,:);
        ddae_rad = radii(pos);
        
        figure(5);
        EV_plot_img(gfp1sl(:,:,track_time+1), swc_ddad); 
        hold on;
        EV_plot_img([], swc_ddae);
        title(num2str(track_time));
        viscircles(centers, radii, 'EdgeColor', 'b');
        viscircles(ddad_cen, ddad_rad, 'EdgeColor', 'r');
        viscircles(ddae_cen, ddae_rad, 'EdgeColor', 'g');
        viscircles(adj_soma(1:2), adj_soma(3), 'EdgeColor', 'c');
        drawnow;
        
        for k = 1:2
            soma_pos{num,k} = zeros(size(gfp1sl,3), 3);
        end
        soma_pos{num,1}(track_time+1,:) = [ddad_cen ddad_rad];
        soma_pos{num,2}(track_time+1,:) = [ddae_cen ddae_rad];
        soma_adj{num} = zeros(size(gfp1sl,3), 3);
        soma_adj{num}(track_time+1,:) = adj_soma;
        
        % track forward
        ddad_cen = soma_pos{num,1}(track_time+1,1:2); 
        ddad_rad = soma_pos{num,1}(track_time+1,3);
        ddae_cen = soma_pos{num,2}(track_time+1,1:2); 
        ddae_rad = soma_pos{num,2}(track_time+1,3);
        soma_adj_cen = adj_soma(1:2);
        soma_adj_rad = adj_soma(3);
        for t = track_time+2:size(gfp1sl,3)
            [centers, radii, metric] = imfindcircles(gfp1sl(:,:,t), soma_opt.radii, ...
                    'ObjectPolarity', 'bright', ...
                    'Sensitivity', soma_opt.sensitiv, ...
                    'EdgeThreshold',soma_opt.edge_thr, ...
                    'Method','TwoStage');
            [ddad_cen, ddad_rad, flag1] = soma_match(gfp1sl(:,:,t), gfp1sl(:,:,t-1), ...
                    centers, radii, ddad_cen, ddad_rad);
            [ddae_cen, ddae_rad, flag2] = soma_match(gfp1sl(:,:,t), gfp1sl(:,:,t-1), ...
                    centers, radii, ddae_cen, ddae_rad);
            [soma_adj_cen, soma_adj_rad, flag3] = soma_match(gfp1sl(:,:,t), gfp1sl(:,:,t-1), ...
                    centers, radii, soma_adj_cen, soma_adj_rad);
            if flag1 || flag2, break; end;
            soma_pos{num,1}(t,:) = [ddad_cen ddad_rad];
            soma_pos{num,2}(t,:) = [ddae_cen ddae_rad];
            soma_adj{num}(t,:) = [soma_adj_cen soma_adj_rad];

            figure(5); imshow(gfp1sl(:,:,t), []); title(t-1);
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(ddad_cen, ddad_rad, 'EdgeColor', 'r');
            viscircles(ddae_cen, ddae_rad, 'EdgeColor', 'g');
            viscircles(soma_adj_cen, soma_adj_rad, 'EdgeColor', 'c');
            drawnow;
        end
        
        % track backward
        ddad_cen = soma_pos{num,1}(track_time+1,1:2); 
        ddad_rad = soma_pos{num,1}(track_time+1,3);
        ddae_cen = soma_pos{num,2}(track_time+1,1:2); 
        ddae_rad = soma_pos{num,2}(track_time+1,3);
        soma_adj_cen = adj_soma(1:2);
        soma_adj_rad = adj_soma(3);
        for t = track_time:-1:1
            [centers, radii, metric] = imfindcircles(gfp1sl(:,:,t), soma_opt.radii, ...
                    'ObjectPolarity', 'bright', ...
                    'Sensitivity', soma_opt.sensitiv, ...
                    'EdgeThreshold',soma_opt.edge_thr, ...
                    'Method','TwoStage');
            [ddad_cen, ddad_rad, flag1] = soma_match(gfp1sl(:,:,t), gfp1sl(:,:,t+1), ...
                    centers, radii, ddad_cen, ddad_rad);
            [ddae_cen, ddae_rad, flag2] = soma_match(gfp1sl(:,:,t), gfp1sl(:,:,t+1), ...
                    centers, radii, ddae_cen, ddae_rad);
            [soma_adj_cen, soma_adj_rad, flag3] = soma_match(gfp1sl(:,:,t), gfp1sl(:,:,t+1), ...
                    centers, radii, soma_adj_cen, soma_adj_rad);
            if flag1 || flag2, break; end;
            soma_pos{num,1}(t,:) = [ddad_cen ddad_rad];
            soma_pos{num,2}(t,:) = [ddae_cen ddae_rad];
            soma_adj{num}(t,:) = [soma_adj_cen soma_adj_rad];

            figure(5); imshow(gfp1sl(:,:,t), []); title(t-1);
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(ddad_cen, ddad_rad, 'EdgeColor', 'r');
            viscircles(ddae_cen, ddae_rad, 'EdgeColor', 'g');
            viscircles(soma_adj_cen, soma_adj_rad, 'EdgeColor', 'c');
            drawnow;
        end
    end
    save([savefile '_somaboth.mat'], 'soma_pos', 'soma_adj');
    disp(['Tracking soma... completed : ' savefile '_somaboth.mat']);
end

% %% use GUI to fit the folding neurons
% results = struct([]);
% for num = 1:length(initswctime)
%     track_time = initswctime(num);
%     swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
%     swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
%     swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
%     swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;
%     
%     [configs,leftpos,rightpos,spline] = ffd_gui_both(gfp, {swc_ddad, swc_ddae}, ...
%             track_time, soma_pos(num,:));
%     %%
%     videoname = sprintf([savefile '_model_both%d.mp4'], num);
%     ffd_plot_both([], gfp, {swc_ddad, swc_ddae}, track_time, soma_pos(num,:), configs, spline);
%     results(num).configs = configs;
%     results(num).leftpos = leftpos;
%     results(num).rightpos = rightpos;
%     results(num).spline = spline;
% end
% save([savefile '_model_both.mat'], 'results');
% disp(['Save tracking to ' savefile '_model_both.mat Done'])

% %% display results
% for num = 1:length(initswctime)
%     track_time = initswctime(num);
%     swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
%     swc(:,3:5) = swc(:,3:5)+1;
%     ffd_plot([], gfp, swc, track_time, soma_pos{num}, results(num).configs, results(num).spline);
% end

% %% model Z-axis
% if exist([savefile '_zflip.mat'], 'file')
%     load([savefile '_zflip.mat'], 'zflip');
% else
%     load([savefile '_model_both.mat'], 'results');
% 
%     for num = 1:length(initswctime)
%         track_time = initswctime(num);
%         swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
%         swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
%         swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
%         swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;
% 
%         zflip = alignZ( gfp, {swc_ddad, swc_ddae}, track_time, soma_pos(num,:), results(num).configs );
%     end
%     save([savefile '_zflip.mat'], 'zflip');
%     disp(['Save tracking to ' savefile '_zflip.mat Done']);
% end

%% display results
load([savefile '_zflip.mat'], 'zflip');
load([savefile '_model_both.mat'], 'results');
load ddaeforward_back.mat
load thetain.mat
opt.tsem = ddaeforward_back_sem;
opt.tmean = ddaeforward_back_mean;
opt.tval = thetain;

for num = 1:length(initswctime)
    opt.zflip = zflip;
    opt.maxdep = 15;
    track_time = initswctime(num);
    swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
    swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
    swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
    swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;

%     videoname = sprintf([savefile '_model_both%d.mp4'], num);
%     ffd_plot_both([], gfp, {swc_ddad, swc_ddae}, track_time, soma_pos(num,:), results(num).configs, results(num).spline, opt);
    
%     videoname = sprintf([savefile '_model_both_stat%d.mp4'], num);
%     theta_ffd_plot_both([], gfp, {swc_ddad, swc_ddae}, track_time, soma_pos(num,:), soma_adj{num}, results(num).configs, results(num).spline, opt);
    videoname = '20170616L7/accordion.mp4';
    ffd_plot_exp(videoname, gfp, {swc_ddad, swc_ddae}, ...
            track_time, {soma_pos{num,1}, soma_pos{num,2}}, ...
            results(num));
end



% %% sanity check
% load([savefile '_model_both.mat'], 'results');
% for num = 1:length(initswctime)
%     numstack = sum(~cellfun(@isempty, results(num).configs(:,1)));
%     Vseq = zeros(size(gfp,1), size(gfp,2), num_slice+1, numstack, 'uint8');
%     track_time = initswctime(num);
%     % read image stack sequence
%     for i = 1:numstack
%         t = track_time + i - 1;
%         for z = 0:num_slice
%             if switchTZ
%                 X = imread( sprintf(gfpfile,z,t) );
%             else
%                 X = imread( sprintf(gfpfile,t,z) );
%             end
%             
%             if normalized
%                 X = single(X);
%                 X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
%             end
%             
%             if adj_slicewise
%                 n = 3;
%                 Idouble = im2double(X);
%                 avg = mean2(Idouble);
%                 sigma = std2(Idouble);
%                 X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
%             end
%             
%             Vseq(:,:,z+1,i) = X;
%         end
%         fprintf('.');
%         if mod(i,100) == 0, fprintf('\n'); end;
%     end
%     
%     swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
%     swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
%     swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
%     swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;
%     
%     %%
%     view3Dmodel(Vseq, track_time, {swc_ddad, swc_ddae}, soma_pos(num,:), results(num).configs, results(num).spline, zflip);
% end



% %% find swc at every timestep
% if exist([savefile '_timestep.mat'], 'file')
%     load([savefile '_timestep.mat'], 'swcs_ddad', 'swcs_ddae');
% else
%     load([savefile '_model_both.mat'], 'results');
% 
%     for num = 1:length(initswctime)
%         track_time = initswctime(num);
%         swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
%         swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
%         swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
%         swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;
% 
%         [swcs_ddad, swcs_ddae] = deformSWC( gfp, {swc_ddad, swc_ddae}, track_time, soma_pos(num,:), results(num).configs, results(num).spline );
%     end
%     save([savefile '_timestep.mat'], 'swcs_ddad', 'swcs_ddae');
%     disp(['Save SWCs to ' savefile '_timestep.mat Done']);
% end
% 
% %% adjust Z
% load([savefile '_model_both.mat'], 'results');
% load([savefile '_timestep.mat'], 'swcs_ddad', 'swcs_ddae');
% for num = 1:length(initswctime)
%     numstack = sum(~cellfun(@isempty, results(num).configs(:,1)));
%     Vseq = zeros(size(gfp,1), size(gfp,2), num_slice+1, numstack, 'uint8');
%     track_time = initswctime(num);
%     % read image stack sequence
%     for i = 1:numstack
%         t = track_time + i - 1;
%         for z = 0:num_slice
%             if switchTZ
%                 X = imread( sprintf(gfpfile,z,t) );
%             else
%                 X = imread( sprintf(gfpfile,t,z) );
%             end
%             
%             if normalized
%                 X = single(X);
%                 X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
%             end
%             
%             if adj_slicewise
%                 n = 3;
%                 Idouble = im2double(X);
%                 avg = mean2(Idouble);
%                 sigma = std2(Idouble);
%                 X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
%             end
%             
%             Vseq(:,:,z+1,i) = X;
%         end
%         fprintf('.');
%         if mod(i,100) == 0, fprintf('\n'); end;
%     end
%     %%
%     close all
%     clearvars himg hddad_e hddad hddae_e hddae
%     col1 = [235 90  235]/255; % ddaD color
%     col2 = [80  200 80 ]/255; % ddaE color
%     [sx,sy,sz] = sphere;
%     rad = 1;%1.5;
%     zgap = 3;
%     [XD, YD, ZD] = meshgrid(1:size(V,2), 1:size(V,1), -2);
%     opt.reg = 2;
%     swc_prev_ddad = swcs_ddad{track_time+1};
%     swc_prev_ddae = swcs_ddae{track_time+1};
%     V_prev = Vseq(:,:,:,1);
%     swc_ddad = cell(1,size(Vseq,4));
%     swc_ddae = cell(1,size(Vseq,4));
%     
%     outputVideo = VideoWriter([savefile '_model_both_adjustZ_zoom.mp4'], 'MPEG-4');
%     open(outputVideo);
%     
%     for ii = 1:size(Vseq,4)
%         t = track_time + ii;
%         if isempty(swcs_ddad{t})
%             break;
%         end
%         swc_ddad{ii} = adjustZ(swcs_ddad{t}, Vseq(:,:,:,ii), swc_prev_ddad, V_prev, opt);
%         swc_ddae{ii} = adjustZ(swcs_ddae{t}, Vseq(:,:,:,ii), swc_prev_ddae, V_prev, opt);
%         
%         if exist('himg', 'var')
%             set(himg, 'CData', gfp(:,:,t));
%         else
%             f = figure(11);
%             set(f, 'Position',[160,200,1050,785],'Color','w');
%             
%             himg = surface(XD,YD,ZD,gfp(:,:,t),...
%                 'FaceColor','texturemap',...
%                 'EdgeColor','none',...
%                 'CDataMapping','direct', ...
%                 'AmbientStrength', 1, ...
%                 'SpecularStrength', 0);
%             colormap(gray(256));
%             hddad = cell(1,size(swc_ddad{ii},1));
%             hddae = cell(1,size(swc_ddae{ii},1));
%             
%             set(gca, 'YDir', 'reverse');
%             lighting gouraud
%             camproj('perspective')
%             axis tight
%             daspect([1 1 1])
%             zlim([-2 10*zgap+2])
%             view([0 50])
%             camlight
%             box on
%         end
%         hold on;
%         % plot ddad
%         nb = swc_ddad{ii}(:,[1 7]);
%         nb(any(nb < 1,2),:) = [];
%         xx = [swc_ddad{ii}(nb(:,1),3), swc_ddad{ii}(nb(:,2),3), nan(size(nb,1),1)]';
%         yy = [swc_ddad{ii}(nb(:,1),4), swc_ddad{ii}(nb(:,2),4), nan(size(nb,1),1)]';
%         zz = [swc_ddad{ii}(nb(:,1),5), swc_ddad{ii}(nb(:,2),5), nan(size(nb,1),1)]';
%         if exist('hddad_e', 'var')
%             set(hddad_e, 'XData', xx(:));
%             set(hddad_e, 'YData', yy(:));
%             set(hddad_e, 'ZData', zgap*zz(:));
%         else
%             hddad_e = plot3(xx(:), yy(:), zgap*zz(:), 'Color', col1, 'LineWidth', 1);
% %             lighting gouraud
% %             camproj('perspective')
% %             axis tight
% %             daspect([1 1 1])
% %             zlim([-2 10*zgap+2])
% %             view([0 25])
% %             camlight
% %             box on
% %             hold on
% %             hddad = cell(1,size(swc_ddad{ii},1));
% %             hddae = cell(1,size(swc_ddae{ii},1));
%         end
%         for i = 1:size(swc_ddad{ii},1)
%             if isempty(hddad{i})
%                 hddad{i} = surf(rad*sx+swc_ddad{ii}(i,3), rad*sy+swc_ddad{ii}(i,4), rad*sz+zgap*swc_ddad{ii}(i,5), ...
%                         'FaceColor', col1, 'EdgeColor', 'none');
%             else
%                 set(hddad{i}, 'XData', rad*sx+swc_ddad{ii}(i,3));
%                 set(hddad{i}, 'YData', rad*sy+swc_ddad{ii}(i,4));
%                 set(hddad{i}, 'ZData', rad*sz+zgap*swc_ddad{ii}(i,5));
%             end
%         end
%         % plot ddae
%         nb = swc_ddae{ii}(:,[1 7]);
%         nb(any(nb < 1,2),:) = [];
%         xx = [swc_ddae{ii}(nb(:,1),3), swc_ddae{ii}(nb(:,2),3), nan(size(nb,1),1)]';
%         yy = [swc_ddae{ii}(nb(:,1),4), swc_ddae{ii}(nb(:,2),4), nan(size(nb,1),1)]';
%         zz = [swc_ddae{ii}(nb(:,1),5), swc_ddae{ii}(nb(:,2),5), nan(size(nb,1),1)]';
%         if exist('hddae_e', 'var')
%             set(hddae_e, 'XData', xx(:));
%             set(hddae_e, 'YData', yy(:));
%             set(hddae_e, 'ZData', zgap*zz(:));
%         else
%             hddae_e = plot3(xx(:), yy(:), zgap*zz(:), 'Color', col2, 'LineWidth', 1);
%         end
%         for i = 1:size(swc_ddae{ii},1)
%             if isempty(hddae{i})
%                 hddae{i} = surf(rad*sx+swc_ddae{ii}(i,3), rad*sy+swc_ddae{ii}(i,4), rad*sz+zgap*swc_ddae{ii}(i,5), ...
%                         'FaceColor', col2, 'EdgeColor', 'none');
%             else
%                 set(hddae{i}, 'XData', rad*sx+swc_ddae{ii}(i,3));
%                 set(hddae{i}, 'YData', rad*sy+swc_ddae{ii}(i,4));
%                 set(hddae{i}, 'ZData', rad*sz+zgap*swc_ddae{ii}(i,5));
%             end
%         end
%         hold off;
%         title(t-1);
%         xlim([-80 80] + swc_ddad{ii}(1,3));
%         writeVideo(outputVideo, getframe(f));
%         pause(0.5);
% %         keyboard;
%         swc_prev_ddad = swc_ddad{ii};
%         swc_prev_ddae = swc_ddae{ii};
%         V_prev = Vseq(:,:,:,ii);
%     end
%     close(outputVideo);
% end
    