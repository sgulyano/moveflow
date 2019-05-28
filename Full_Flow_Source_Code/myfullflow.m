functionname='myfullflow.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
currentdir = pwd;
cd(functiondir);

addpath(genpath('external'));%load external libraries

if ~exist('debug', 'dir'), mkdir('debug'); end
if ~exist('../results/', 'dir'), mkdir('../results/'); end

close all

%%
% FullFlow parameters
ratio=3;%downsample ratio
maxDisp=100;%smaller displacement for fast processing
%maxDisp=242;%maximal displacement used in the evaluation
opt=struct('setting','sintel','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',3);%Sintel paramters
%opt=struct('setting','kitti','outOfRange',0.22,'lambda',0.021,'truncation',1e8,'maxIter',3,'inverse_b',0);%Kitti paramters

% Our parameters
dataset = 1;
switch dataset
    case 1
        directory = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/';
        filename = 'larva3_2_z%d_t%03d.tif';
        savefile = '../larva3_frames/larva3_fullflow';
        start_stack = 0;
        num_stacks = 8;
        num_timestep = 409;
        adj_thr = [0.02,0.25];
    case 2
        directory = '~/Desktop/LipingsData/1121larva1_1/';
        filename = '1122larva1_1good_z%d_t%03d.tif';
        savefile = '../results/1122larva1_fullflow';
        start_stack = 0;
        num_stacks = 8;
        num_timestep = 751;
        adj_thr = [0.0,0.1];
    case 3
        directory = '~/Desktop/LipingsData/1202DLa2_1good_subset/';
        filename = '1202DLa2_1good_Subset_z%d_t%03d.tif';
        savefile = '../results/1202DLa2_fullflow';
        start_stack = 3;
        num_stacks = 7;
        num_timestep = 899;
        adj_thr = [0.17,0.7];
end
vector_space = 8;
SAVE_VIDEO = true;

%% find optical flow
info = imfinfo(sprintf(fullfile(directory, filename), 0, 0));

U = zeros(info.Height,info.Width,num_timestep);
V = zeros(info.Height,info.Width,num_timestep);
Is = zeros(info.Height,info.Width,num_timestep+1);

vecmask = false(size(U,1), size(U,2));
vecmask(vector_space:vector_space:size(U,1),vector_space:vector_space:size(U,2)) = true;

im1 = [];
h1 = [];
for t = 0:num_timestep
    disp(['Time step: ' num2str(t)]);
    % read image stack
    Xs = [];
    for z = start_stack:num_stacks
        fname = sprintf(fullfile(directory, filename), z, t);
        X = imread(fname);
        X = single(X);
        X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        Xs(:,:,z+1) = X;
    end
    I = max(Xs,[],3);
    Iadj = imadjust(uint8(I),adj_thr,[0,1]);
    Iadj = double(Iadj) ./ 255;
    Is(:,:,t+1) = Iadj;
%     figure(1); imshow(Iadj,[]); drawnow;
    im2 = repmat(Iadj, 1, 1, 3);
    if ~isempty(im1)
        % apply Full Flow
        forward=fullflow(im1,im2,ratio,maxDisp,opt);%compute the forward flow
        backward=fullflow(im2,im1,ratio,maxDisp,opt);%compute the backward flow
        removeOcclusion(forward,backward,ratio,'debug/match.txt');%save the matches to match.txt after forward-backward consistency checking
        flow=runEpicflow(im1,im2,'debug/match.txt',opt.setting);%apply Epicflow interpolation on the matches
        
        Is(:,:,t) = Iadj;
        U(:,:,t) = flow(:,:,1);
        V(:,:,t) = flow(:,:,2);
        
        sc = max(max(sqrt(U(:,:,t).^2+V(:,:,t).^2)));
        figure(2);
        if isempty(h1)
            [h1, h2] = plotFlow(U(:,:,t), V(:,:,t), Is(:,:,t), vector_space, sc/2);
        else
            set(h1, 'CData', Is(:,:,t));
            Utmp = U(:,:,t); Utmp(~vecmask) = 0;
            set(h2, 'UData', Utmp);
            Vtmp = V(:,:,t); Vtmp(~vecmask) = 0;
            set(h2, 'VData', Vtmp);
            set(h2, 'AutoScaleFactor', sc/2);
        end
        title(num2str(t));
        drawnow;
    end
    im1 = im2;
end

save([savefile '.mat'], 'U', 'V', 'Is', '-v7.3');

%% plot result and record in video
if SAVE_VIDEO,
    outputVideo = VideoWriter([savefile '.avi']);
    open(outputVideo);
    
    vecmask = false(size(U,1), size(U,2));
    vecmask(vector_space:vector_space:size(U,1),vector_space:vector_space:size(U,2)) = true;
    h1 = [];
    for t = 1:num_timestep
%         figure(1), clf, subplot(2,2,1); imshow(Is(:,:,t-1));
%         subplot(2,2,2); imshow(Is(:,:,t));
%         subplot(2,2,3); imshow(imfuse(Is(:,:,t-1),Is(:,:,t)));
%         subplot(2,2,4); imshow(flowToColor(cat(3, U(:,:,t), V(:,:,t))));
        
        sc = max(max(sqrt(U(:,:,t).^2+V(:,:,t).^2)));
        figure(2);
        if isempty(h1)
            [h1, h2] = plotFlow(U(:,:,t), V(:,:,t), Is(:,:,t), vector_space, sc/2);
        else
            set(h1, 'CData', Is(:,:,t));
            Utmp = U(:,:,t); Utmp(~vecmask) = 0;
            set(h2, 'UData', Utmp);
            Vtmp = V(:,:,t); Vtmp(~vecmask) = 0;
            set(h2, 'VData', Vtmp);
            set(h2, 'AutoScaleFactor', sc/2);
        end
        title(num2str(t));
        
        set(gcf, 'Position', [0 0 800 400]);
        writeVideo(outputVideo, getframe(gcf));
    end
    close(outputVideo);
end

cd(currentdir);