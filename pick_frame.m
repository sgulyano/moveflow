%%
dataset = 1;
SAVE_VIDEO = true;
switch dataset
    case 1
        directory = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/';
        filename = 'larva3_2_z%d_t%03d.tif';
        savefile = 'larva3';
        num_stacks = 8;
        num_timestep = 409;
        adj_thr = [0.0,0.08];
        frames = [198, 201, 229, 252, 340, 342];
    case 2
        directory = '~/Desktop/LipingsData/1121larva1_1/';
        filename = '1122larva1_1good_z%d_t%03d.tif';
        savefile = '1122larva1';
        num_stacks = 8;
        num_timestep = 751;
        adj_thr = [0.0,0.1];
        frames = [43, 47, 184, 285, 445, 555, 660];
    case 3
        directory = '~/Desktop/LipingsData/1202DLa2_1good_subset/';
        filename = '1202DLa2_1good_Subset_z%d_t%03d.tif';
        savefile = '1202DLa2';
        num_stacks = 7;
        num_timestep = 899;
        adj_thr = [0.17,0.7];
end
%%
disp('Read Image:');
V = [];
for t = 1:num_timestep
    Xs = [];
    for z = 0:num_stacks
        fname = sprintf(fullfile(directory, filename), z, t);
        X = imread(fname);
        X = imadjust(X,adj_thr,[0,1]);
        X = double(X);
        X = X-min(X(:)); X = X/max(X(:));
        Xs(:,:,z+1) = X;
    end
    V(:,:,t) = max(Xs,[],3)*255;
    fprintf('.');
    if mod(t,100) == 0
        fprintf('\n');
    end
end
fprintf('\n');
figure(1), imshow3D(V);

%%
if ~exist([savefile '_frames'], 'dir')
    mkdir([savefile '_frames']);
end
for frame = frames
    frame_dir = [savefile '_frames/' savefile '_t' num2str(frame)];
    if ~exist(frame_dir, 'dir')
        mkdir(frame_dir);
    end
    Xs = [];
    for z = 0:num_stacks
        fname = sprintf(fullfile(directory, filename), z, frame);
        X = imread(fname);
        X = imadjust(X,adj_thr,[0,1]);
        X = double(X);
        X = X-min(X(:)); X = X/max(X(:));
        Xs(:,:,z+1) = X;
        sname = sprintf(fullfile(frame_dir, filename), z, frame);
        imwrite(X, sname);
    end
    figure(2); clf; imshow3D(Xs);
%     keyboard;
end