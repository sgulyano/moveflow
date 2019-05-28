close all; clear all; clc;

addpath(genpath('Full_Flow_Source_Code'))
addpath(genpath('toolbox'))
addpath(genpath('UGM'))

dataset = 1;
opt.DEBUG = false;
SAVE_VIDEO = true;
switch dataset
    case 1
        filename = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/larva3_2_z%d_t%03d.tif';
        savefile = 'larva3_frames/larva3';
        savename = '~/Desktop/snakegraphmodel/swc_dan/ZstackL1_3_2grayscale/larva3_2_t%03d.swc';
        num_timestep = 409;
        num_stacks = 8;
        imadj_param = [0 0.25];
        Isize = [192 512 num_stacks+1];
    case 2
        directory = '~/Desktop/LipingsData/1121larva1_1/';
        filename = '1122larva1_1good_z%d_t%03d.tif';
        savefile = '1122larva1';
        savename = '~/Desktop/snakegraphmodel/swc_dan/1122Larva1/1122Larva1_t%03d.swc';
        num_timestep = 751;
        num_stacks = 8;
        imadj_param = [0 0.1];
        Isize = [192 512 num_stacks+1];
    case 3
        directory = '~/Desktop/LipingsData/1202DLa2_1good_subset/';
        filename = '1202DLa2_1good_Subset_z%d_t%03d.tif';
        savefile = '1202DLa2';
        savename = '~/Desktop/snakegraphmodel/swc_dan/1202DLa2/1202DLa2_t%03d.swc';
        num_timestep = 899;
        num_stacks = 7;
end

%% load initial trace
% read vector field
load([savefile '_fullflow.mat'])
load([savefile '_soma.mat'], 'soma_pos');

%% track neuron by moving the trace using MRF
% read trace
swclist = dir([savefile '_t*.swc']);
initswctime = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
initswctime = cellfun(@(x)str2double(x{1}), initswctime);
% start tracking
traces = cell(num_timestep+1, length(initswctime));
for num = 1:length(initswctime)
    track_time = initswctime(num);
    disp(['Tracking : ' savefile '_t' num2str(track_time) '.swc']);
    swc = read_swc_file( [savefile '_t' num2str(track_time) '.swc'] );
    swc(:,3:5) = swc(:,3:5)+1; swc(:,5) = min(max(swc(:,5), 1), num_stacks+1);

    % track forward
    traces{track_time+1,num} = swc;
    for t = track_time+1:num_timestep
        swc = traces{t,num};
        tic;
        [swc2, flag] = mrftrack(swc, Is(:,:,t+1), soma_pos{num}(t+1,:), soma_pos{num}(t,:), ...
                traces{track_time+1,num}, soma_pos{num}(track_time+1,:), opt);
        toc;
        if flag, break; end
        traces{t+1,num} = swc2;
        figure(4); EV_plot_img( Is(:,:,t+1), traces{t+1,num} ); 
        title(num2str(t)); drawnow;
    end
    % track backward
    for t = track_time:-1:1
        swc = traces{t+1,num};
        tic;
        [swc2, flag] = mrftrack(swc, Is(:,:,t), soma_pos{num}(t,:), soma_pos{num}(t+1,:), ...
                traces{track_time+1,num}, soma_pos{num}(track_time+1,:), opt);
        toc;
        if flag, break; end
        traces{t,num} = swc2;
        figure(4); EV_plot_img( Is(:,:,t), traces{t,num} ); title(num2str(t-1)); drawnow;
    end
end
save([savefile '_track.mat'], 'traces', 'initswctime');
disp(['Tracking neuron... completed : ' savefile '_track.mat']);

%% save tracking result to video
if SAVE_VIDEO
    img = zeros(size(Is,1),size(Is,2), num_timestep+1, 'uint8');
    for t = 0:num_timestep
        Xs = uint8([]);
        for z = 0:num_stacks
            Xs(:,:,z+1) = imread( sprintf(filename, z, t) );
        end
        img(:,:,t+1) = max(Xs,[],3);
    end
    [ a_sum, a_mu, a_num ] = track2vdo( [savefile '_track.avi'], img, traces );
    track2csv( [savefile '_track.csv'], soma_pos, num_timestep, a_sum, a_mu, a_num );
end
