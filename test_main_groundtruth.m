addpath(genpath('toolbox/utils'))
addpath(genpath('fix_slice'))

dataset = 1;

switch dataset
    case 1
        directory = 'C:\Users\y_sar\Documents\MATLAB\DanData\1110DLa_2_P1';
        filename = '1110DLa_2_P1_z%d_t%02d.tif';
        savedir = 'C:\Users\y_sar\Documents\moveflow\1110DLa_2_P1';
end

load(fullfile(savedir, 'user_input.mat'), ...
        'params', 'soma_cen', 'soma_rad', 'forwardmovement', ...
        'initswctime', 'ddad_swc', 'ddae_swc', 'dd_plane');

num_pairs = length(ddad_swc);
            
% read the middle frame for tuning parameter
[num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);

[gfp, gfp1sl] = readImg(directory, filename, params);
           

%% track a pair of ddaD and ddaE by fitting a plane through folding
for i = 1:num_pairs
    dd_plane(i).configs = [];
end

for i = 1:num_pairs
    [configs, zflip, spline] = ffd_gui_both_groundtruth(gfp, {ddad_swc{i}, ddae_swc{i}}, ...
            initswctime(i), {soma_cen{1,i}, soma_cen{2,i}}, dd_plane(i).leftpos, dd_plane(i).rightpos);
    dd_plane(i).configs = configs;
    dd_plane(i).spline = spline;
    dd_plane(i).zflip = zflip;
end

save(fullfile(savedir, 'groundtruth.mat'), ...
        'params', 'soma_cen', 'soma_rad', 'forwardmovement', ...
        'initswctime', 'ddad_swc', 'ddae_swc', 'dd_plane');

for i = 1:num_pairs
    opt.zflip = dd_plane(i).zflip;
    plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)], opt);
end

disp('Done');