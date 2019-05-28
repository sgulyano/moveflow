dataset = 1;
SAVE_VIDEO = true;

Options.SigmaFluid = 4;
Options.Verbose = 0;

switch dataset
    case 1
        directory = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/';
        filename = 'larva3_2_z%d_t%03d.tif';
        savefile = 'larva3';
        start_stack = 5;
        num_stacks = 8;
        num_timestep = 409;
        adj_thr = [0.0,0.1];
        mag_thr = 1000;
    case 2
        directory = '~/Desktop/LipingsData/1121larva1_1/';
        filename = '1122larva1_1good_z%d_t%03d.tif';
        savefile = '1122larva1';
        start_stack = 5;
        num_stacks = 8;
        num_timestep = 751;
        adj_thr = [0.0,0.05];
        mag_thr = 1000;
    case 3
        directory = '~/Desktop/LipingsData/1202DLa2_1good_subset/';
        filename = '1202DLa2_1good_Subset_z%d_t%03d.tif';
        savefile = '1202DLa2_fullflow';
        start_stack = 3;
        num_stacks = 7;
        num_timestep = 899;
        adj_thr = [0.17,0.7];
end

%% select frame to correct
load(['results/' savefile '_fullflow.mat']);
SqrtMag = sqrt(squeeze(sum(sum(U.^2+V.^2,1),2)));

figure(3); plot(SqrtMag);
hold on; plot([1 num_timestep], [1000 1000], 'r'); hold off;

t_idx = SqrtMag > mag_thr;
t_idx(t_idx(2:end)) = 1;
t_idx([false; t_idx(1:end-1)]) = 1;
disp(['Total number of stacks to fix alignment = ' num2str(sum(t_idx))]);
clearvars SqrtMag U V Is

%% correct alignment
t_idx_arr = find(t_idx);
Vreg_tmp = cell(1,length(t_idx_arr));

st_time = tic;
poolobj = parpool(4);
parfor n = 1:length(t_idx_arr)
    t = t_idx_arr(n);
    disp(['Timestep: ', num2str(t)]);
    % read image stack
    V = [];
    for z = 0:num_stacks
        fname = sprintf(fullfile(directory, filename), z, t);
        X = imread(fname);
        X = double(X);
        X = X-min(X(:)); X = X/max(X(:));
        V(:,:,z+1) = X;
    end
    
    % adjust intensity for registration
    Vadj = V;
    for z = 1:size(V,3)
        Vadj(:,:,z) = imadjust(V(:,:,z),adj_thr,[0,1]);
    end
    figure(2), subplot(2,1,1); imshow(max(Vadj,[],3));

    % register along depth
    for i = start_stack:num_stacks;
        [Vadj(:,:,i+1), Bx, By] = register_images(Vadj(:,:,i+1),Vadj(:,:,i),Options);
        V(:,:,i+1) = movepixels(V(:,:,i+1),Bx,By,[], 3);
    end
    for i = start_stack:-1:2;
        [Vadj(:,:,i-1), Bx, By] = register_images(Vadj(:,:,i-1),Vadj(:,:,i),Options);
        V(:,:,i-1) = movepixels(V(:,:,i-1),Bx,By,[], 3);
    end
    
    Vreg_tmp{n} = uint8(V*255);

    subplot(2,1,2); imshow(max(Vadj,[],3)); drawnow;
end
delete(poolobj);
disp('Total time');
toc(st_time);

Vregis = cell(1, num_timestep);
Vregis(t_idx_arr) = Vreg_tmp;
save(['../results/' savefile '_regis.mat'], 'Vregis', '-v7.3');

%% save to video
if SAVE_VIDEO,
    load(['../results/' savefile '_regis.mat']);    
    
    outputVideo = VideoWriter(['../results/' savefile '_regis.avi']);
    open(outputVideo);
    
    for t = 0:num_timestep
        % load original image stack
        V = uint8([]);
        for z = 0:num_stacks
            fname = sprintf(fullfile(directory, filename), z, t);
            X = imread(fname);
            X = double(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(X*255);
            V(:,:,z+1) = X;
        end
        % aligned image stack
        if t > 0 && ~isempty(Vregis{t})
            Vreg = Vregis{t};
        else
            Vreg = V;
        end
        % adjust intensity for display
        Vadj = V;
        Vreg_adj = Vreg;
        for z = 0:num_stacks
            Vadj(:,:,z+1) = imadjust(V(:,:,z+1),[0 0.25],[0,1]);
            Vreg_adj(:,:,z+1) = imadjust(Vreg(:,:,z+1),[0 0.25],[0,1]);
        end
        % plot
        figure(1); subplot(2,1,1); imshow(max(Vadj,[],3));
        title(num2str(t));
        subplot(2,1,2); imshow(max(Vreg_adj,[],3));
        writeVideo(outputVideo, getframe(gcf));
    end
    close(outputVideo);
end
