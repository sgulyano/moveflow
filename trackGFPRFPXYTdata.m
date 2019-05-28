close all; clear all; clc;

addpath(genpath('Full_Flow_Source_Code'))
addpath(genpath('toolbox'));
addpath(genpath('UGM'));

functionname = 'trackGFPRFPXYTdata.m';
functiondir = which(functionname);
functiondir = functiondir(1:end-length(functionname));

%% user parameters
SAVE_VIDEO = true;
dataset = 1;
switch dataset
    case 1
        gfpfile = '~/Desktop/LipingsData/GFPRFPXYTdata/s7 tiff green/Larva 4_Series007_Crop001_Crop001_t%04d_ch00.tif';
        rfpfile = '~/Desktop/LipingsData/GFPRFPXYTdata/s7 tiff red/Larva 4_Series007_Crop001_Crop002_t%04d_ch00.tif';
        num_img = 815;
        savefile = 'Larva4s007/Larva4s007';
        gfpadjrange = [0.04 0.5];
        rfpadjrange = [0 .1];
        pickframes = {240:245, 645:660, 790:810};
        sigcur = 5;
        soma_opt.adjval = [.2 .5];
        soma_opt.radii = [7 15];
        soma_opt.sensitiv = .9;
        soma_opt.edge_thr = 0.1;
        soma_opt.scale = 1;
end

%% read data
disp('Reading GFP and RFP images');
gfpinfo = imfinfo(sprintf(gfpfile,0));
rfpinfo = imfinfo(sprintf(rfpfile,0));
if gfpinfo.Height ~= rfpinfo.Height || gfpinfo.Width ~= rfpinfo.Width,
    error('GFP and RFP size mismatch');
end
gfp = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
rfp = zeros(rfpinfo.Height, rfpinfo.Width, num_img+1, 'uint8');
for num = 0:num_img
    gfp(:,:,num+1) = imread( sprintf(gfpfile,num) );
    rfp(:,:,num+1) = imread( sprintf(rfpfile,num) );
    fprintf('.');
    if mod(num+1,100)==0, fprintf('\n'); end
end
fprintf('\n');

if exist([savefile '_fullflow.mat'], 'file')
    load([savefile '_fullflow.mat']);
else
    %% preprocessing
    disp('Preprocessing images...');
    gfpmed = medfilt3(gfp);
    rfpmed = medfilt3(rfp);
    gfpadj = gfpmed;
    rfpadj = rfpmed;
    for num = 0:num_img
        gfpadj(:,:,num+1) = imadjust(gfpmed(:,:,num+1), gfpadjrange, [0 1]);
        rfpcur = NfdctDemo(double(rfpmed(:,:,num+1)), sigcur);
        rfpadj(:,:,num+1) = imadjust(max(rfpcur,0) ./ 255, rfpadjrange, [0 1]) * 255;
    end
    figure(1); imshow3D(gfpadj);
%     figure(2); imshow3D(rfpadj);

    %% pick frames for manual tracing (using neuTube)
    i = 0;  % save one image as stack at a time where i = 1:length(pickframes) 
            % so we can read in neuTube
    if i > 0
        frame_i = max(rfpadj(:,:,pickframes{i}),[],3);
        imwrite(frame_i, [savefile '_frame1.tif']);
        imwrite(frame_i, [savefile '_frame2.tif']);
    end
    
    %% compute Optical Flow (using FullFlow)
    disp('Compute Optical Flow using FullFlow...');
    cd([functiondir '/Full_Flow_Source_Code']);
    addpath(genpath('external'));       %load external libraries
    if ~exist('debug', 'dir'), mkdir('debug'); end

    % FullFlow parameters
    ratio=3;%downsample ratio
    maxDisp=100;%smaller displacement for fast processing
    %maxDisp=242;%maximal displacement used in the evaluation
    opt=struct('setting','sintel','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',3,'parfor',true);%Sintel paramters
    %opt=struct('setting','kitti','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',0);%Kitti paramters

    U = zeros(size(gfpadj)-[0 0 1]);
    V = zeros(size(gfpadj)-[0 0 1]);

    vecmask = false(size(U,1), size(U,2));
    vector_space = 8;
    vecmask(vector_space:vector_space:size(U,1),vector_space:vector_space:size(U,2)) = true;

    flowst = tic;
    parfor t = 1:num_img
        disp(['Time step: ' num2str(t)]);
        im1 = repmat(im2double(gfpadj(:,:,t)), [1 1 3]);
        im2 = repmat(im2double(gfpadj(:,:,t+1)), [1 1 3]);

        % apply Full Flow
        forward=fullflow(im1,im2,ratio,maxDisp,opt);%compute the forward flow
        backward=fullflow(im2,im1,ratio,maxDisp,opt);%compute the backward flow
        removeOcclusion(forward,backward,ratio,'debug/match.txt');%save the matches to match.txt after forward-backward consistency checking
        flow=runEpicflow(im1,im2,'debug/match.txt',opt.setting,t);%apply Epicflow interpolation on the matches

        U(:,:,t) = flow(:,:,1);
        V(:,:,t) = flow(:,:,2);

        if ~opt.parfor
            sc = max(max(sqrt(U(:,:,t).^2+V(:,:,t).^2)));
            figure(3);
            if ~exist('h1', 'var') || ~ishandle(h1)
                [h1, h2] = plotFlow(U(:,:,t), V(:,:,t), im2, vector_space, sc/2);
            else
                set(h1, 'CData', im2);
                Utmp = U(:,:,t); Utmp(~vecmask) = 0;
                set(h2, 'UData', Utmp);
                Vtmp = V(:,:,t); Vtmp(~vecmask) = 0;
                set(h2, 'VData', Vtmp);
                set(h2, 'AutoScaleFactor', sc/2);
            end
            title(num2str(t));
            drawnow;
        end
    end
    stat.flow = toc(flowst);
    if opt.parfor
        delete debug/*.bin debug/*.png debug/*.flo
    end
    cd(functiondir);

    save([savefile '_fullflow.mat'], 'U', 'V', 'gfpadj', 'rfpadj', '-v7.3');
    disp(['Compute vector flow... completed : ' savefile '_fullflow.mat']);
    
    %% write to video
    if SAVE_VIDEO
        imgflow2vdo( [savefile '_gfp.avi'], gfpadj );
        imgflow2vdo( [savefile '_rfp.avi'], rfpadj );
        imgflow2vdo( [savefile '_fullflow.avi'], gfpadj, U, V );
    end
end


%% get swc list
swclist = dir([savefile '_t*.swc']);
initswctime = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
initswctime = cellfun(@(x)str2double(x{1}), initswctime);

if exist([savefile '_soma.mat'], 'file')
    load([savefile '_soma.mat']);
else
    %% track soma using hough transform
    gfpmed = medfilt3(gfp);
    figure(4); clf; imshow3D(gfpmed);

    %% track soma for each trace
    disp('Tracking soma... started');
    [xx,yy]=meshgrid(1:size(gfp,2),1:size(gfp,1));
    soma = cell(length(initswctime),1);
    soma_pos = cell(length(initswctime),1);
    for num = 1:length(initswctime)        
        disp(['Track : ' swclist(num).name]);
        track_time = initswctime(num);
        swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
        swc(:,3:5) = swc(:,3:5)+1;

        % find initial soma location
        [ centers, radii ] = detect_soma( gfpmed(:,:,track_time+1), soma_opt );
        [~, pos] = min(sum(bsxfun(@minus, centers, swc(1,3:4)).^2,2));
        soma_cen = centers(pos,:);
        soma_rad = radii(pos);

        figure(5);  EV_plot_img(gfpadj(:,:,track_time+1), swc);
        viscircles(centers, radii, 'EdgeColor', 'b');
        viscircles(soma_cen, soma_rad, 'EdgeColor', 'r');
        drawnow;
        
        soma_pos{num} = zeros(size(gfpmed,3), 3);
        soma_pos{num}(track_time+1,:) = [soma_cen soma_rad];
        soma{num} = false(size(gfpmed));
        soma{num}(:,:,track_time+1) = ((xx-soma_cen(1)).^2+(yy-soma_cen(2)).^2) <= soma_rad^2;
        % track forward
        centroid = soma_cen;
        for t = track_time+1:num_img
            [centers, radii, ~] = detect_soma( gfpmed(:,:,t+1), soma_opt );
            [centroid, rad, flag] = soma_hough(U(:,:,t), V(:,:,t), centers, radii, centroid);
            if flag, break; end;
            soma{num}(:,:,t+1) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
            soma_pos{num}(t+1,:) = [centroid rad];

            figure(5); imshow( gfpadj(:,:,t+1), [] ); title(num2str(t)); drawnow;
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(centroid, rad, 'EdgeColor', 'r');
            drawnow;
        end
        % track backward
        centroid = soma_cen;
        for t = track_time:-1:1
            [centers, radii, ~] = detect_soma( gfpmed(:,:,t), soma_opt );
            [centroid, rad, flag] = soma_hough(-U(:,:,t), -V(:,:,t), centers, radii, centroid);
            if flag, break; end;
            soma{num}(:,:,t) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
            soma_pos{num}(t,:) = [centroid rad];

            figure(5); imshow( gfpadj(:,:,t), [] ); title(num2str(t-1)); drawnow;
            viscircles(centers, radii, 'EdgeColor', 'b');
            viscircles(centroid, rad, 'EdgeColor', 'r');
            drawnow;
        end
    end
    
    save([savefile '_soma.mat'], 'soma', 'soma_pos');
    disp(['Tracking soma... completed : ' savefile '_soma.mat']);
end

%% track neuron by moving the trace using MRF
if exist([savefile '_track.mat'], 'file')
    load([savefile '_track.mat']);
else
    if isa(gfpadj,'uint8')
        gfpadj = im2double(gfpadj);
    end
    traces = cell(num_img+1, length(initswctime));
    for num = 1:length(initswctime)
        track_time = initswctime(num);
        disp(['Tracking : ' savefile '_t' num2str(track_time) '.swc']);
        swc = read_swc_file( sprintf('%s_t%03d.swc', savefile, track_time) );
        swc(:,3:5) = swc(:,3:5)+1; swc(:,5) = min(max(swc(:,5), 1), 2);

        % track forward
        traces{track_time+1,num} = swc;
        for t = track_time+1:num_img
            swc = traces{t,num};
            [swc2, flag] = mrftrack(swc, gfpadj(:,:,t+1), soma_pos{num}(t+1,:), soma_pos{num}(t,:), ...
                    traces{track_time+1,num}, soma_pos{num}(track_time+1,:));
            if flag, break; end
            traces{t+1,num} = swc2;
            figure(6); EV_plot_img( gfpadj(:,:,t+1), traces{t+1,num} ); 
            title(num2str(t)); drawnow;
        end
        % track backward
        for t = track_time:-1:1
            swc = traces{t+1,num};
            [swc2, flag] = mrftrack(swc, gfpadj(:,:,t), soma_pos{num}(t,:), soma_pos{num}(t+1,:), ...
                        traces{track_time+1,num}, soma_pos{num}(track_time+1,:));
            if flag, break; end
            traces{t,num} = swc2;
            figure(6); EV_plot_img( gfpadj(:,:,t), traces{t,num} ); title(num2str(t-1)); drawnow;
        end
    end
    save([savefile '_track.mat'], 'traces', 'initswctime');
    disp(['Tracking neuron... completed : ' savefile '_track.mat']);
    
    %% save tracking result to video
    if SAVE_VIDEO
        track_opt.soma_pos = soma_pos;
        [a_sum, a_mu, a_num] = track2vdo( [savefile '_track.avi'], gfp, traces, track_opt );
        track2csv( [savefile '_track.csv'], soma_pos, num_img, a_sum, a_mu, a_num );
    end
end
