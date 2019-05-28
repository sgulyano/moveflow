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
dataset = 6;
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
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        soma_opt = struct();                            % options for soma detector
        soma_opt.radii = [15 30];
        soma_opt.sensitiv = 0.6;
        soma_opt.edge_thr = 0.3;
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

if exist([savefile '_fullflow.mat'], 'file')
    load([savefile '_fullflow.mat']);
else
    %% preprocessing
    disp('Preprocessing images...');
    if strcmp(use_medfilt,'3D')
        gfpmed = medfilt3(gfp);
    elseif strcmp(use_medfilt,'2D')
        gfpmed = gfp;
        for num = 0:num_img
            gfpmed(:,:,num+1) = medfilt2(gfp(:,:,num+1));
        end
    else
        gfpmed = gfp;
    end
    gfpadj = gfpmed;
    for num = 0:num_img
        gfpadj(:,:,num+1) = imadjust(gfpmed(:,:,num+1), gfpadjrange, [0 1]);
    end
    
    if DEBUG, figure(1); imshow3D(gfpadj); end
    if ~isempty(rfpfile)
        rfpmed = medfilt3(rfp);
        rfpadj = rfpmed;
        for num = 0:num_img
            rfpadj(:,:,num+1) = imadjust(rfpmed(:,:,num+1), rfpadjrange, [0 1]);
        end
        if DEBUG, figure(2); imshow3D(rfpadj); end
    else
        rfpadj = [];
    end

%     %% pick frames for manual tracing (using neuTube)
%     i = 1;  % save one image as stack at a time where i = 1:length(pickframes) 
%             % so we can read in neuTube
%     if i > 0
%         t = pickframes(i);
%         for z = 0:num_slice
%             frame_i = imread( sprintf(gfpfile,z,t) );
%             X = frame_i;
%             n = 2;
%             Idouble = im2double(X);
%             avg = mean2(Idouble);
%             sigma = std2(Idouble);
%             X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
%             frame_i = X;
%             frame_i = imadjust(frame_i, gfpadjrange, [0 1]);
%             imwrite(frame_i, [savefile '_frame' num2str(z) '.tif']);
%         end
%     end
% 
%     %% compute Optical Flow (using FullFlow)
%     disp('Compute Optical Flow using FullFlow...');
%     cd([functiondir '/Full_Flow_Source_Code']);
%     addpath(genpath('external'));       %load external libraries
%     if ~exist('debug', 'dir'), mkdir('debug'); end
% 
%     % FullFlow parameters
%     ratio=3;%downsample ratio
%     maxDisp=100;%smaller displacement for fast processing
%     %maxDisp=242;%maximal displacement used in the evaluation
%     opt=struct('setting','sintel','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',3,'parfor',false);%Sintel paramters
%     %opt=struct('setting','kitti','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',0);%Kitti paramters
% 
%     U = zeros(size(gfpadj)-[0 0 1]);
%     V = zeros(size(gfpadj)-[0 0 1]);
% 
%     vecmask = false(size(U,1), size(U,2));
%     vector_space = 8;
%     vecmask(vector_space:vector_space:size(U,1),vector_space:vector_space:size(U,2)) = true;
% 
%     flowst = tic;
%     for t = 1:num_img
%         disp(['Time step: ' num2str(t)]);
%         im1 = repmat(im2double(gfpadj(:,:,t)), [1 1 3]);
%         im2 = repmat(im2double(gfpadj(:,:,t+1)), [1 1 3]);
% 
%         % apply Full Flow
%         forward=fullflow(im1,im2,ratio,maxDisp,opt);%compute the forward flow
%         backward=fullflow(im2,im1,ratio,maxDisp,opt);%compute the backward flow
%         removeOcclusion(forward,backward,ratio,'debug/match.txt');%save the matches to match.txt after forward-backward consistency checking
%         flow=runEpicflow(im1,im2,'debug/match.txt',opt.setting,t);%apply Epicflow interpolation on the matches
% 
%         U(:,:,t) = flow(:,:,1);
%         V(:,:,t) = flow(:,:,2);
% 
%         if ~opt.parfor && DEBUG
%             sc = max(max(sqrt(U(:,:,t).^2+V(:,:,t).^2)));
%             figure(3);
%             if ~exist('h1', 'var') || ~ishandle(h1)
%                 [h1, h2] = plotFlow(U(:,:,t), V(:,:,t), im2, vector_space, sc/2);
%             else
%                 set(h1, 'CData', im2);
%                 Utmp = U(:,:,t); Utmp(~vecmask) = 0;
%                 set(h2, 'UData', Utmp);
%                 Vtmp = V(:,:,t); Vtmp(~vecmask) = 0;
%                 set(h2, 'VData', Vtmp);
%                 set(h2, 'AutoScaleFactor', sc/2);
%             end
%             title(num2str(t));
%             drawnow;
%         end
%     end
%     stat.flow = toc(flowst);
%     if opt.parfor
%         delete debug/*.bin debug/*.png debug/*.flo
%     end
%     cd(functiondir);
%     
%     save([savefile '_fullflow.mat'], 'U', 'V', 'gfpadj', 'rfpadj', '-v7.3');
%     disp(['Compute vector flow... completed : ' savefile '_fullflow.mat']);
%     
%     %% write to video
%     if SAVE_VIDEO
%         % plot image processing result
%         imgflow2vdo( [savefile '_gfp.avi'], gfpadj );
%     %     imgflow2vdo( [savefile '_rfp.avi'], rfpadj );
%         imgflow2vdo( [savefile '_fullflow.avi'], gfpadj, U, V );
%     end
end

%% get swc list
swclist = dir([savefile '_t*.swc']);
initswctime = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
initswctime = cellfun(@(x)str2double(x{1}), initswctime);

if exist([savefile '_soma.mat'], 'file')
    load([savefile '_soma.mat']);
else
    %% track soma using hough transform
    if strcmp(use_medfilt,'3D')
        gfp1slmed = medfilt3(gfp1sl);
    elseif strcmp(use_medfilt,'2D')
        gfp1slmed = gfp1sl;
        for num = 0:num_img
            gfp1slmed(:,:,num+1) = medfilt2(gfp1sl(:,:,num+1));
        end
    else
        gfp1slmed = gfp1sl;
    end
    
    if DEBUG, figure(4); imshow3D(gfp1slmed); end

    % track soma for each trace
    disp('Tracking soma... started');
    [xx,yy]=meshgrid(1:size(gfp1sl,2),1:size(gfp1sl,1));
    soma = cell(length(initswctime),1);
    soma_pos = cell(length(initswctime),1);
    for num = 1:length(initswctime)
        disp(['Track : ' swclist(num).name]);
        track_time = initswctime(num);
        swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
        swc(:,3:5) = swc(:,3:5)+1;

        soma_opt.radii = [20 40]
        soma_opt.sensitiv = 0.7;
        soma_opt.edge_thr = 0.3;
        % find initial soma location
        [ centers, radii ] = detect_soma( gfp1slmed(:,:,track_time+1), soma_opt );
        [~, pos] = min(sum(bsxfun(@minus, centers, swc(1,3:4)).^2,2));
        soma_cen = centers(pos,:);
        soma_rad = radii(pos);

        figure(5);  EV_plot_img(gfpadj(:,:,track_time+1), swc); title(num2str(track_time));
        viscircles(centers, radii, 'EdgeColor', 'b');
        viscircles(soma_cen, soma_rad, 'EdgeColor', 'r');
        drawnow;
        
        soma_pos{num} = zeros(size(gfp1slmed,3), 3);
        soma_pos{num}(track_time+1,:) = [soma_cen soma_rad];
        soma{num} = false(size(gfp1slmed));
        soma{num}(:,:,track_time+1) = ((xx-soma_cen(1)).^2+(yy-soma_cen(2)).^2) <= soma_rad^2;
        % track forward
        centroid = soma_cen; rad = soma_rad;
        for t = track_time+1:num_img
            [centers, radii, ~] = detect_soma( gfp1slmed(:,:,t+1), soma_opt );
            if use_normxcorr
                [centroid, rad, flag] = soma_match(gfp1slmed(:,:,t+1), gfp1slmed(:,:,t), centers, radii, centroid, rad);
            else
                [centroid, rad, flag] = soma_hough(U(:,:,t), V(:,:,t), centers, radii, centroid);
            end
            if flag, break; end;
            soma{num}(:,:,t+1) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
            soma_pos{num}(t+1,:) = [centroid rad];

            figure(5); imshow( gfpadj(:,:,t+1), [] ); title(num2str(t));
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(centroid, rad, 'EdgeColor', 'r');
            drawnow;
        end
        % track backward
        centroid = soma_cen; rad = soma_rad;
        for t = track_time:-1:1
            [centers, radii, ~] = detect_soma( gfp1slmed(:,:,t), soma_opt );
            
            if use_normxcorr
                [centroid, rad, flag] = soma_match(gfp1slmed(:,:,t), gfp1slmed(:,:,t+1), centers, radii, centroid, rad);
            else
                [centroid, rad, flag] = soma_hough(-U(:,:,t), -V(:,:,t), centers, radii, centroid);
            end
            if flag, break; end;
            soma{num}(:,:,t) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
            soma_pos{num}(t,:) = [centroid rad];

            figure(5); imshow( gfpadj(:,:,t), [] ); title(num2str(t-1));
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(centroid, rad, 'EdgeColor', 'r');
            drawnow;
        end
    end
    save([savefile '_soma.mat'], 'soma', 'soma_pos');
    disp(['Tracking soma... completed : ' savefile '_soma.mat']);
end

%% track neuron by moving the trace using MRF
if FULLFLOW,
    trackfile = [savefile '_track_fullflow'];
else
    trackfile = [savefile '_track'];
end
if exist([trackfile '.mat'], 'file')
    load([trackfile '.mat']);
else
    if isa(gfpadj,'uint8')
        gfpadj = im2double(gfpadj);
    end
    traces = cell(num_img+1, length(initswctime));
    for num = 1:length(initswctime)
        track_time = initswctime(num);
        disp(['Tracking : ' savefile '_t' num2str(track_time) '.swc']);
        swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
        swc(:,3:5) = swc(:,3:5)+1; swc(:,5) = min(max(swc(:,5), 1), num_slice+1);

        % track forward
        traces{track_time+1,num} = swc;
        for t = track_time+1:num_img
            swc = traces{t,num};
            if FULLFLOW,
                [swc2, flag] = opticalflowtrack( swc, U(:,:,t), V(:,:,t) );
            else
                tic;
                [swc2, flag] = mrftrack(swc, gfpadj(:,:,t+1), soma_pos{num}(t+1,:), soma_pos{num}(t,:), ...
                        traces{track_time+1,num}, soma_pos{num}(track_time+1,:), track_opt);
                toc;
            end
            if flag, break; end
            traces{t+1,num} = swc2;
            figure(6); EV_plot_img( gfpadj(:,:,t+1), traces{t+1,num} ); 
            title(num2str(t)); drawnow;
        end
        % track backward
        for t = track_time:-1:1
            swc = traces{t+1,num};
            if FULLFLOW,
                [swc2, flag] = opticalflowtrack( swc, -U(:,:,t), -V(:,:,t) );
            else
                tic;
                [swc2, flag] = mrftrack(swc, gfpadj(:,:,t), soma_pos{num}(t,:), soma_pos{num}(t+1,:), ...
                            traces{track_time+1,num}, soma_pos{num}(track_time+1,:), track_opt);
                toc;
            end
            if flag, break; end
            traces{t,num} = swc2;
            figure(6); EV_plot_img( gfpadj(:,:,t), traces{t,num} ); title(num2str(t-1)); drawnow;
        end
    end
    save([trackfile '.mat'], 'traces', 'initswctime');
    disp(['Tracking neuron... completed : ' trackfile '.mat']);
    
    %% save tracking result to video
    if SAVE_VIDEO
        track_opt.soma_pos = soma_pos;
        [a_sum, a_mu, a_num] = track2vdo( [trackfile '.avi'], gfp, traces, track_opt );
        track2csv( [trackfile '.csv'], soma_pos, num_img, a_sum, a_mu, a_num );
    end
end
