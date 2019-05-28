addpath(genpath('toolbox/utils'))
% load('20180404-L6-B/user_input.mat');
% [gfp, gfp1sl] = readImg(directory, filename, params);
% %%
% for i = 1:length(dd_plane)
%     plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)] );
% end


load('20180404-L7-B/user_input.mat');
directory = '~/Desktop/LipingsData/20180404mcdGFP_Curvature/20180404-L7-B';
[gfp, gfp1sl] = readImg(directory, filename, params);
%%
for i = 1:length(dd_plane)
    plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)] );
end

%%
ddad_lefts = zeros(1,num_pairs);
for i = 1:num_pairs
    ddad_lefts(i) = sign(ddad_swc{i}(1,3) - ddae_swc{i}(1,3));
end
ddaD_Left = mode(ddad_lefts) == -1;

% convert to excel            
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
[~, directions] = cellfun(@(x)linegrowing(x(:,1)), soma_cen, 'UniformOutput', false);
opt.directions = directions;
convert2excel('20180404-L7-B', gfp, size(gfp,3)-1, initswctime, ddad_swc, ddae_swc, soma_cen, soma_rad, dd_plane, opt);
