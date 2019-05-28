close all; clear all; clc;

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
        gfpfile = '~/Desktop/LipingsData/20170616L7_align/20170616L7_align_z%d_t%03d.tif';     % Ca images
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
        adj_soma = [29, 86];
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

% %% read img data
% disp('Reading GFP and RFP images');
% Vinfo = imfinfo(sprintf(gfpfile,0,0));
% 
% gfp = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, motion_end_frames-pickframes+1, 'uint8');
% V = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, motion_end_frames-pickframes+1, 'uint8');
% for t = pickframes:motion_end_frames
%     for z = 0:num_slice
%         if switchTZ
%             X = imread( sprintf(gfpfile,z,t) );
%         else
%             X = imread( sprintf(gfpfile,t,z) );
%         end
%         if normalized
%             X = single(X);
%             X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
%         end
%         gfp(:,:,z+1,t-pickframes+1) = X;
%         
%         if adj_slicewise
%             n = 3;
%             Idouble = im2double(X);
%             avg = mean2(Idouble);
%             sigma = std2(Idouble);
%             X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
%         end
%         V(:,:,z+1,t-pickframes+1) = X;
%     end
% end
% ZNUM = start_slice+1:end_slice+1;
% gfp = gfp(:,:,ZNUM,:);
% V = V(:,:,ZNUM,:);
% fprintf('\n');
% figure(1), imshow3D(reshape(V,[128,512,size(V,3)*size(V,4)]))

%% get swc list
swclist = dir([savefile '_t*.swc']);
initswctime = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
initswctime = cellfun(@(x)str2double(x{1}), initswctime);

load([savefile '_soma.mat']);
% if exist([savefile '_soma_slice.mat'], 'file')
%     load([savefile '_soma_slice.mat']);
% else
%     %% track soma using hough transform
%     % track soma for each trace
%     disp('Tracking soma... started');
%     soma_pos = cell(length(initswctime),1);
%     for num = 1%:length(initswctime)
%         track_time = initswctime(num);
%         swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
%         swc(:,3:5) = swc(:,3:5)+1;
% 
%         [centers, radii, metric] = imfindcircles(gfp(:,:,pickslice,1), soma_opt.radii, ...
%                 'ObjectPolarity', 'bright', ...
%                 'Sensitivity', soma_opt.sensitiv, ...
%                 'EdgeThreshold',soma_opt.edge_thr, ...
%                 'Method','TwoStage');
%         
%         [~, pos] = min(sum(bsxfun(@minus, centers, swc(1,3:4)).^2,2));
%         soma_cen = centers(pos,:);
%         soma_rad = radii(pos);
% 
%         figure(5);  EV_plot_img(V(:,:,pickslice,1), swc); title(num2str(track_time));
%         viscircles(centers, radii, 'EdgeColor', 'b');
%         viscircles(soma_cen, soma_rad, 'EdgeColor', 'r');
%         drawnow;
%         
%         soma_pos{num} = zeros(size(V,4), size(V,3), 3);
%         soma_pos{num}(1,:,:) = repmat([soma_cen soma_rad], size(V,3), 1);
%         
%         % track forward
%         centroid = repmat(soma_cen, size(V,3), 1); rad = repmat(soma_rad, size(V,3), 1);
%         for t = 2:size(V,4)
%             for z = 1:size(V,3)
%                 [centers, radii, metric] = imfindcircles(gfp(:,:,z,t), soma_opt.radii, ...
%                         'ObjectPolarity', 'bright', ...
%                         'Sensitivity', soma_opt.sensitiv, ...
%                         'EdgeThreshold',soma_opt.edge_thr, ...
%                         'Method','TwoStage');
%                 [centroid(z,:), rad(z), flag] = soma_match(gfp(:,:,z,t), gfp(:,:,z,t-1), ...
%                         centers, radii, centroid(z,:), rad(z));
%                 if flag, break; end;
%                 soma_pos{num}(t,z,:) = [centroid(z,:) rad(z)];
%                 
%                 figure(5); imshow( V(:,:,z,t), [] ); title(num2str([t+track_time-1,z]));
%                 viscircles(centers, radii, 'EdgeColor', 'b');
%                 viscircles(centroid(z,:), rad(z), 'EdgeColor', 'r');
%                 drawnow;
%             end
%         end
%         
%     end
%     save([savefile '_soma_slice.mat'], 'soma_pos');
%     disp(['Tracking soma... completed : ' savefile '_soma_slice.mat']);
% end


% %% use GUI to fit the folding neurons
% if exist([savefile '_model.mat'], 'file')
%     load([savefile '_model.mat'], 'results');
% else
%     results = struct([]);
%     for num = 1:length(initswctime)
%         track_time = initswctime(num);
%         swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
%         swc(:,3:5) = swc(:,3:5)+1;
%         [configs,leftpos,rightpos,spline] = ffd_gui(gfp, swc, ...
%                 track_time, soma_pos{num});
% 
%         videoname = sprintf([savefile '_model%d.avi'], num);
%         ffd_plot(videoname, gfp, swc, track_time, soma_pos{num}, configs, spline);
%         results(num).configs = configs;
%         results(num).leftpos = leftpos;
%         results(num).rightpos = rightpos;
%         results(num).spline = spline;
%     end
%     save([savefile '_model.mat'], 'results');
%     disp(['Save tracking to ' savefile '_model.mat Done'])
% end 

% %% display results
% for num = 1:length(initswctime)
%     track_time = initswctime(num);
%     swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
%     swc(:,3:5) = swc(:,3:5)+1;
%     ffd_plot([], gfp, swc, track_time, soma_pos{num}, results(num).configs, results(num).spline);
% end

%% plot X-axis FFD radius and delta F from Gcamp data
load([savefile '_model.mat'], 'results');
load thetain.mat
opt.tval = thetain;

switch('ddaDbw')
    case 'ddaDbw'
        load ddadbackward_back.mat
        opt.tsem = ddadbackward_back_sem;
        opt.tmean = ddadbackward_back_mean;
    case 'ddaEfw'
        load ddaeforward_back.mat
        opt.tsem = ddaeforward_back_sem;
        opt.tmean = ddaeforward_back_mean;
    otherwise
        opt.tsem = [];
        opt.tmean = [];
end


for num = 1:length(initswctime)
    track_time = initswctime(num);
    swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
    swc(:,3:5) = swc(:,3:5)+1;
    
    videoname = sprintf([savefile '_model%d.mp4'], num);
    ffd_plot(videoname, gfp, swc, track_time, soma_pos{num}, results(num).configs, results(num).spline);
    videoname = sprintf([savefile '_model_stat%d.mp4'], num);
    theta_ffd_plot(videoname, gfp, swc, track_time, soma_pos{num}, soma_adj{num}, results(num).configs, results(num).spline, opt);
end