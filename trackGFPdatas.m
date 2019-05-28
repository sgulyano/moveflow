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
        pickframes = [48, 140, 332, 458, 633];          % frame picked for manual tracing
        pickslice = 5;                                  % image slice for tracking soma
        use_medfilt = [];                               % use median filter in soma detection
        
        soma_opt = [struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1),...
                struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1),...
                struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1),...
                struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1),...
                struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1)];
        adj_soma = [161, 67, 6; ...
                    140, 71, 5; ...
                    130, 80, 5; ...
                    436, 71, 5; ...
                    479, 63, 6];
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

% %% pick frames for manual tracing (using neuTube)
% for t = pickframes
%     dirname = sprintf('%s_t%03d', savefile, t);
%     mkdir( dirname );
%     for z = 0:num_slice
%         frame_i = imread( sprintf(gfpfile,z,t) );
%         X = frame_i;
%         n = 2;
%         Idouble = im2double(X);
%         avg = mean2(Idouble);
%         sigma = std2(Idouble);
%         X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
%         frame_i = X;
%         frame_i = imadjust(frame_i, gfpadjrange, [0 1]);
%         imwrite(frame_i, [dirname '/frame' num2str(z) '.tif']);
%     end
% end

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
    
%     if DEBUG, figure(4); imshow3D(gfp1slmed); end

    % track soma for each trace
    disp('Tracking soma... started');
    [xx,yy]=meshgrid(1:size(gfp1sl,2),1:size(gfp1sl,1));
    soma = cell(length(initswctime),1);
    soma_pos = cell(length(initswctime),1);
    soma_adj = cell(length(initswctime),1);
    for num = 1:length(initswctime)
        disp(['Track : ' swclist(num).name]);
        track_time = initswctime(num);
        swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
        swc(:,3:5) = swc(:,3:5)+1;

%         soma_opt(num).radii = [5 8];
%         soma_opt(num).sensitiv = 0.99;
%         soma_opt(num).edge_thr = 0.05;
        % find initial soma location
        [centers, radii, ~] = imfindcircles(gfp1slmed(:,:,track_time+1), soma_opt(num).radii, ...
                'ObjectPolarity', 'bright', ...
                'Sensitivity', soma_opt(num).sensitiv, ...
                'EdgeThreshold', soma_opt(num).edge_thr, ...
                'Method','TwoStage');
        [~, pos] = min(sum(bsxfun(@minus, centers, swc(1,3:4)).^2,2));
        soma_cen = centers(pos,:);
        soma_rad = radii(pos);
        
        %[~, pos] = min(sum(bsxfun(@minus, centers, adj_soma).^2,2));
        soma2_cen = adj_soma(num,1:2);%centers(pos,:);
        soma2_rad = adj_soma(num,3);%radii(pos);
        

        figure(5);  EV_plot_img(gfp1sl(:,:,track_time+1), swc); title(num2str(track_time));
        viscircles(centers, radii, 'EdgeColor', 'b');
        viscircles(soma_cen, soma_rad, 'EdgeColor', 'r');
        viscircles(soma2_cen, soma2_rad, 'EdgeColor', 'g');
        drawnow;
%         keyboard;
        soma_pos{num} = zeros(size(gfp1slmed,3), 3);
        soma_pos{num}(track_time+1,:) = [soma_cen soma_rad];
        soma{num} = false(size(gfp1slmed));
        soma{num}(:,:,track_time+1) = ((xx-soma_cen(1)).^2+(yy-soma_cen(2)).^2) <= soma_rad^2;
        soma_adj{num} = zeros(size(gfp1slmed,3), 3);
        soma_adj{num}(track_time+1,:) = [soma2_cen soma2_rad];

        % track forward
        centroid = soma_cen; rad = soma_rad;
        centroid2 = soma2_cen; rad2 = soma2_rad;
        for t = track_time+1:num_img
            [centers, radii, ~] = imfindcircles(gfp1slmed(:,:,t+1), soma_opt(num).radii, ...
                    'ObjectPolarity', 'bright', ...
                    'Sensitivity', soma_opt(num).sensitiv, ...
                    'EdgeThreshold', soma_opt(num).edge_thr, ...
                    'Method','TwoStage');
            if use_normxcorr
                [centroid, rad, flag] = soma_match(gfp1slmed(:,:,t+1), gfp1slmed(:,:,t), centers, radii, centroid, rad);
                [centroid2, rad2, f2] = soma_match(gfp1slmed(:,:,t+1), gfp1slmed(:,:,t), centers, radii, centroid2, rad2);
                if f2
                    centroid2 = [0 0]; rad2 = 0;
                end
            else
                [centroid, rad, flag] = soma_hough(U(:,:,t), V(:,:,t), centers, radii, centroid);
            end
%             if t > 100, keyboard; end
            if flag, break; end;
            soma{num}(:,:,t+1) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
            soma_pos{num}(t+1,:) = [centroid rad];
            soma_adj{num}(t+1,:) = [centroid2 rad2];

            figure(5); imshow( gfp1sl(:,:,t+1), [] ); title(num2str(t));
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(centroid, rad, 'EdgeColor', 'r');
            viscircles(centroid2, rad2, 'EdgeColor', 'g');
            drawnow;
        end
        % track backward
        centroid = soma_cen; rad = soma_rad;
        centroid2 = soma2_cen; rad2 = soma2_rad;
        for t = track_time:-1:1
            [centers, radii, ~] = imfindcircles(gfp1slmed(:,:,t), soma_opt(num).radii, ...
                    'ObjectPolarity', 'bright', ...
                    'Sensitivity', soma_opt(num).sensitiv, ...
                    'EdgeThreshold', soma_opt(num).edge_thr, ...
                    'Method','TwoStage');
            if use_normxcorr
                [centroid, rad, flag] = soma_match(gfp1slmed(:,:,t), gfp1slmed(:,:,t+1), centers, radii, centroid, rad);
                [centroid2, rad2, f2] = soma_match(gfp1slmed(:,:,t), gfp1slmed(:,:,t+1), centers, radii, centroid2, rad2);
                if f2
                    centroid2 = [0 0]; rad2 = 0;
                end
            else
                [centroid, rad, flag] = soma_hough(-U(:,:,t), -V(:,:,t), centers, radii, centroid);
            end
            if flag, break; end;
            soma{num}(:,:,t) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
            soma_pos{num}(t,:) = [centroid rad];
            soma_adj{num}(t,:) = [centroid2 rad2];

            figure(5); imshow( gfp1sl(:,:,t), [] ); title(num2str(t-1));
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(centroid, rad, 'EdgeColor', 'r');
            viscircles(centroid2, rad2, 'EdgeColor', 'g');
            drawnow;
%             keyboard;
        end
    end
    save([savefile '_soma.mat'], 'soma', 'soma_pos', 'soma_adj');
    disp(['Tracking soma... completed : ' savefile '_soma.mat']);
%     keyboard;
end

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

        baseline_dist = abs(soma_pos{num}(track_time+1,1) - soma_adj{num}(track_time+1,1));
        % track forward
        traces{track_time+1,num} = swc;
        for t = track_time+1:num_img
            if soma_adj{num}(t+1,1) == 0
                soma_dist = baseline_dist;
            else
                soma_dist = abs(soma_pos{num}(t+1,1) - soma_adj{num}(t+1,1));
            end
            swc = traces{t,num};
            if FULLFLOW,
                [swc2, flag] = opticalflowtrack( swc, U(:,:,t), V(:,:,t) );
            else
                tic;
                model_trace = traces{track_time+1,num};
                xfix = model_trace(1,3);
                ratio = soma_dist / baseline_dist;
                model_trace(:,3) = ((model_trace(:,3)-xfix) * ratio)+ xfix;
                [swc2, flag] = mrftrack(swc, gfpadj(:,:,t+1), soma_pos{num}(t+1,:), soma_pos{num}(t,:), ...
                        model_trace, soma_pos{num}(track_time+1,:), track_opt);
                toc;
            end
            if flag, break; end
            traces{t+1,num} = swc2;
            figure(6); EV_plot_img( gfpadj(:,:,t+1), traces{t+1,num} ); 
            title(num2str(t)); drawnow;
        end
        % track backward
        for t = track_time:-1:1
            soma_dist = abs(soma_pos{num}(t,1) - soma_adj{num}(t,1));
            swc = traces{t+1,num};
            if FULLFLOW,
                [swc2, flag] = opticalflowtrack( swc, -U(:,:,t), -V(:,:,t) );
            else
                tic;
                model_trace = traces{track_time+1,num};
                xfix = model_trace(1,3);
                ratio = soma_dist / baseline_dist;
                model_trace(:,3) = ((model_trace(:,3)-xfix) * ratio)+ xfix;
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
        track_opt.soma_adj = soma_adj;
        [a_len, neu_len, len_baseline] = track2vdo_length( [trackfile '_len.avi'], gfp, traces, initswctime, track_opt );
        
        T = table((0:num_img)');
        T.Properties.VariableNames{end} = 'ImageNumber';
        Ts = cell(1,length(soma_pos));
        Tadj = cell(1,length(soma_pos));
        for i = 1:length(soma_pos)
            Ts{i} = table(soma_pos{i}(:,1), soma_pos{i}(:,2));
            Ts{i}.Properties.VariableNames = {['Xsoma_' num2str(i)], ['Ysoma_' num2str(i)]};
            Tadj{i} = table(soma_adj{i}(:,1), soma_adj{i}(:,2));
            Tadj{i}.Properties.VariableNames = {['Xprevsoma_' num2str(i)], ['Yprevsoma_' num2str(i)]};
        end
        Tlen = table(a_len);
        Tlen.Properties.VariableNames{1} = 'dD_D';
        Tneu = table(neu_len);
        Tneu.Properties.VariableNames{1} = 'TotalDendriteLength';
        Tflat = table(len_baseline);
        Tflat.Properties.VariableNames{1} = 'FlatDendriteLength';
        writetable([T Ts{:} Tadj{:} Tlen Tneu Tflat], [trackfile '_len.csv']);
        disp(['Save track result to ' trackfile '_len.csv' ' ... Done']);
    end
end
