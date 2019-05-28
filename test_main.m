addpath(genpath('toolbox/utils'))
addpath(genpath('fix_slice'))

savedir = fullfile(pwd, '20180404-L6-B');
load(fullfile(savedir, 'user_input.mat'))

directory = 'C:\Users\y_sar\Documents\MATLAB\DanData\20180404mcdGFP_Curvature\20180404-L6-B_aligned';
filename = '20180404-L6-B_z%02d_t%02d.tif';

% tuning parameter
[num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);
[gfp, gfp1sl] = readImg(directory, filename, params);
%%
num_pairs =  length(initswctime);
for i = 1:num_pairs
    plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), 'Test' );
end

% %% track a pair of ddaD and ddaE by fitting a plane through folding
% num_pairs =  length(initswctime);
% results = struct([]);
% for i = 1:num_pairs
%     [configs, leftpos, rightpos] = fit_plane_gui_both(gfp, {ddad_swc{i}, ddae_swc{i}}, ...
%             initswctime(i), {soma_cen{1,i}, soma_cen{2,i}});
% %     videoname = fullfile(savedir, sprintf('model_both%d.mp4', i));
% %     ffd_plot_both(videoname, gfp, {ddad_swc{i}, ddae_swc{i}}, ...
% %             initswctime(i), {soma_cen{1,i}, soma_cen{2,i}}, ...
% %             configs, spline, struct('is_gui', true));
%     results(i).configs = configs;
%     results(i).leftpos = leftpos;
%     results(i).rightpos = rightpos;
%     results(i).spline = spline;
%     
%     plot_track( gfp, results(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)] );
%    
% end

%% convert to excel            
A = load('thetain.mat', 'thetain');
opt.is_gui = true;
opt.tval = A.thetain;
opt.ddad_left = ddaD_Left;
if forwardmovement
    A = load('ddaeforward_back.mat', 'ddaeforward_back_sem', 'ddaeforward_back_mean');
    opt.tsem = A.ddaeforward_back_sem;
    opt.tmean = A.ddaeforward_back_mean;
else
    A = load('ddadbackward_back.mat', 'ddadbackward_back_sem', 'ddadbackward_back_mean');
    opt.tsem = A.ddadbackward_back_sem;
    opt.tmean = A.ddadbackward_back_mean;
end

% [ lines, direction ] = linegrowing( soma_cen{1,2}(:,1), struct('debug', true) );

[~, directions] = cellfun(@(x)linegrowing(x(:,1)), soma_cen, 'UniformOutput', false);
opt.directions = directions;
convert2excel([], gfp, num_img, initswctime, ddad_swc, ddae_swc, soma_cen, soma_rad, dd_plane, opt);
