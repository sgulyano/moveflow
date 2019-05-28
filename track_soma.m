close all; clear all; clc;

addpath(genpath('Full_Flow_Source_Code'))
addpath(genpath('demon_registration_version_8f'));
addpath('toolbox/AMT');
% addpath(genpath('AOSLevelsetSegmentationToolboxM'));
addpath('toolbox/utils')

dataset = 1;
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
        filename = '~/Desktop/LipingsData/1121larva1_1/1122larva1_1good_z%d_t%03d.tif';
        savefile = '1122larva1_frames/1122larva1';
        savename = '~/Desktop/snakegraphmodel/swc_dan/1122Larva1/1122Larva1_t%03d.swc';
        num_timestep = 751;
        num_stacks = 8;
        imadj_param = [0 0.1];
        Isize = [192 512 num_stacks+1];
    case 3
        filename = '~/Desktop/LipingsData/1202DLa2_1good_subset/1202DLa2_1good_Subset_z%d_t%03d.tif';
        savefile = '1202DLa2/1202DLa2';
        savename = '~/Desktop/snakegraphmodel/swc_dan/1202DLa2/1202DLa2_t%03d.swc';
        num_timestep = 899;
        num_stacks = 7;
end

%% load initial trace
% read vector field
load([savefile '_fullflow.mat'])
% read trace
flist = dir([savefile '_t*.swc']);

timestep = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), flist);
timestep = cellfun(@(x)str2double(x{1}), timestep);

%% read image
info = imfinfo( sprintf(filename, 0, 0) );
imgs = zeros(info.Height,info.Width,num_timestep+1);

disp('Read Images');
for t = 0:num_timestep
    % read image stack
    Xs = [];
    for z = 4%0:num_stacks
        fname = sprintf(filename, z, t);
        X = imread(fname);
        X = single(X);
        X = X-min(X(:)); X = X/max(X(:));
        Xs(:,:,z+1) = X;
    end
    imgs(:,:,t+1) = Xs(:,:,5);%max(Xs,[],3);
    fprintf('.');
    if mod(t+1,100)==0, fprintf('\n'); end
end
fprintf('\n');

%% track soma using hough transform
[xx,yy]=meshgrid(1:size(imgs,2),1:size(imgs,1));
soma = cell(length(flist),1);
soma_pos = cell(length(flist),1);
for num = 1:length(flist)
    disp(['Track : ' flist(num).name]);
    track_time = timestep(num);
    swc = read_swc_file( [savefile '_t' num2str(track_time) '.swc'] );
    swc(:,3:5) = swc(:,3:5)+1;
    
    % find initial soma location
    [ centers, radii ] = detect_soma( imgs(:,:,track_time+1) );
    [~, pos] = min(sum(bsxfun(@minus, centers, swc(1,3:4)).^2,2));
    soma_cen = centers(pos,:);
    soma_rad = radii(pos);
    
    figure(2); subplot(3,3,num); imshow( Is(:,:,track_time+1) );
    viscircles(centers, radii, 'EdgeColor', 'b');
    viscircles(soma_cen, soma_rad, 'EdgeColor', 'r');
    hold on; scatter(swc(1,3), swc(1,4), 100, 'r+');
    drawnow;
    
    
    soma_pos{num} = zeros(size(Is,3), 3);
    soma_pos{num}(track_time+1,:) = [soma_cen soma_rad];
    soma{num} = false(size(imgs));
    soma{num}(:,:,track_time+1) = ((xx-soma_cen(1)).^2+(yy-soma_cen(2)).^2) <= soma_rad^2;
    % track forward
    centroid = soma_cen;
    for t = track_time+1:num_timestep
        [centers, radii, ~] = detect_soma( imgs(:,:,t+1) );
        [centroid, rad, flag] = soma_hough(U(:,:,t), V(:,:,t), centers, radii, centroid);
        if flag, break; end;
        soma{num}(:,:,t+1) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
        soma_pos{num}(t+1,:) = [centroid rad];
        
        figure(3); imshow( Is(:,:,t+1), [] ); title(num2str(t)); drawnow;
        hold on;
        plot(centroid(1), centroid(2), 'rx');
        contour(soma{num}(:,:,t+1), [.5 .5], 'r');
        hold off;
        drawnow;
    end
    % track backward
    centroid = soma_cen;
    for t = track_time:-1:1
        [centers, radii, ~] = detect_soma( imgs(:,:,t) );
        [centroid2, rad, flag] = soma_hough(-U(:,:,t), -V(:,:,t), centers, radii, centroid);
        if flag, break; end;
%         if t < 187
%             keyboard;
%         end
        centroid = centroid2;
        soma{num}(:,:,t) = ((xx-centroid(1)).^2+(yy-centroid(2)).^2) <= rad^2;
        soma_pos{num}(t,:) = [centroid rad];
        
        figure(3); imshow( Is(:,:,t), [] ); title(num2str(t-1)); drawnow;
        hold on;
        plot(centroid(1), centroid(2), 'rx');
        contour(soma{num}(:,:,t), [.5 .5], 'r');
        hold off;
        drawnow;
    end
end
save([savefile '_soma.mat'], 'soma', 'soma_pos');

%% track soma using snake
% soma = cell(length(flist),1);
% for num = 1:length(flist)
%     disp(['Track : ' flist(num).name]);
%     track_time = timestep(num);
%     
%     soma{num} = false(size(U) + [0 0 1]);
%     soma{num}(:,:,track_time+1) = poly2mask(init_soma(num).position(:,1), init_soma(num).position(:,2), size(U,1), size(U,2));
%     % track forward
%     centroid = init_soma(num).centroid;
%     snake = AC_remesh_close(init_soma(num).position,1,'equal','close');
%     for t = track_time+1:num_timestep
%         [snake, centroid, flag] = soma_snake(snake, U(:,:,t), V(:,:,t), imgs(:,:,t+1), centroid, 1);
%         soma{num}(:,:,t) = poly2mask(snake(:,1), snake(:,2), size(U,1), size(U,2));
%         if flag, break; end;
%         
%         figure(2); imshow( imgs(:,:,t+1), [] ); title(num2str(t)); drawnow;
%         hold on;
%         plot(centroid(1), centroid(2), 'rx');
%         plot(snake(:,1), snake(:,2), 'r-');
%         hold off;
%         drawnow;
%     end
%     % track backward
%     centroid = init_soma(num).centroid;
%     snake = AC_remesh_close(init_soma(num).position,1,'equal','close');    
%     for t = track_time:-1:1
%         [snake, centroid, flag] = soma_snake(snake, U(:,:,t), V(:,:,t), imgs(:,:,t), centroid, -1);
%         soma{num}(:,:,t) = poly2mask(snake(:,1), snake(:,2), size(U,1), size(U,2));
%         if flag, break; end;
%         
%         figure(2); imshow( imgs(:,:,t), [] ); title(num2str(t-1)); drawnow;
%         hold on;
%         plot(centroid(1), centroid(2), 'rx');
%         plot(snake(:,1), snake(:,2), 'r-o');
%         hold off;
%         drawnow;
%     end
% end

%% segment soma in each frame using level set
% soma = cell(length(flist),1);
% for num = 1:length(flist)
%     track_time = timestep(num);
%     
%     soma{num} = false(size(U) + [0 0 1]);
%     soma{num}(:,:,track_time+1) = poly2mask(init_soma(num).position(:,1), init_soma(num).position(:,2), size(U,1), size(U,2));
%     % track forward
%     BW = soma{num}(:,:,track_time+1);
%     centroid = init_soma(num).centroid;
%     for t = track_time+1:num_timestep
%         [ centroid, BW, flag ] = soma_levelset( U(:,:,t), V(:,:,t), imgs(:,:,t+1), centroid, BW, 1 );
%         soma{num}(:,:,t) = BW;
%         if flag, break; end;
%         
%         figure(2); imshow( imgs(:,:,t+1), [] ); title(num2str(t)); drawnow;
%         hold on;
%         plot(centroid(1), centroid(2), 'rx');
%         contour(BW, [.5 .5], 'r');
%         hold off;
%         drawnow;
%     end
%     % track backward
%     BW = soma{num}(:,:,track_time+1);
%     centroid = init_soma(num).centroid;
%     for t = track_time:-1:1
%         [ centroid, BW, flag ] = soma_levelset( U(:,:,t), V(:,:,t), imgs(:,:,t), centroid, BW, -1 );
%         soma{num}(:,:,t) = BW;
%         if flag, break; end;
%                 
%         figure(2); imshow( imgs(:,:,t), [] ); title(num2str(t-1)); drawnow;
%         hold on;
%         plot(centroid(1), centroid(2), 'rx');
%         contour(BW, [.5 .5], 'r');
%         hold off;
%         drawnow;
%         pause;
%     end
% end

%% track soma using only optical flow
% [xx, yy] = meshgrid(1:size(U,2), 1:size(U,1));
% soma = cell(length(flist),1);
% for num = 1:length(flist)
%     track_time = timestep(num);
%     
%     soma{num} = false(size(U)+[0 0 1]);
%     soma{num}(:,:,track_time+1) = poly2mask(init_soma(num).position(:,1), init_soma(num).position(:,2), size(U,1), size(U,2));
%     % track forward
%     BW = soma{num}(:,:,track_time+1);
%     centroid = init_soma(num).centroid;
%     for t = track_time+1:num_timestep
%         dx = interp2(U(:,:,t), centroid(1), centroid(2));
%         dy = interp2(V(:,:,t), centroid(1), centroid(2));
%         if isnan(dx) || isnan(dy), break; end
%         centroid = centroid + [dx dy];
%         
%         [ii, jj] = ind2sub(size(BW), find(BW));
%         nx = round(jj + dx); ny = round(ii + dy);
%         idxOut = ny < 1 | ny > size(U,1) | nx < 1 | nx > size(U,2);
%         
%         idx = sub2ind(size(BW), ny(~idxOut), nx(~idxOut));
%         BW = false(size(BW));
%         BW(idx) = true;
%         
%         soma{num}(:,:,t+1) = BW;
%         
%         if sum(idxOut) > length(idxOut)/2
%             break
%         end
%         
%         figure(2); imshow( Is(:,:,t+1) ); title(num2str(t)); drawnow;
%         hold on;
%         plot(centroid(1), centroid(2), 'rx');
%         contour(BW, [.5 .5], 'r');
%         hold off;
%         drawnow;
%     end
%     % track backward
%     BW = soma{num}(:,:,track_time+1);
%     centroid = init_soma(num).centroid;
%     for t = track_time:-1:1
%         % since U,V are forward mapping, do interpolate scattered data
%         nx = xx + U(:,:,t);
%         ny = yy + V(:,:,t);
%         
%         dx = griddata(nx, ny, -U(:,:,t), centroid(1), centroid(2));
%         dy = griddata(nx, ny, -V(:,:,t), centroid(1), centroid(2));
%         if isnan(dx) || isnan(dy), break; end
%         centroid = centroid + [dx dy];
%         
%         [ii, jj] = ind2sub(size(BW), find(BW));
%         nx = round(jj + dx); ny = round(ii + dy);
%         idxOut = ny < 1 | ny > size(U,1) | nx < 1 | nx > size(U,2);
%         
%         idx = sub2ind(size(BW), ny(~idxOut), nx(~idxOut));
%         BW = false(size(BW));
%         BW(idx) = true;
%         
%         soma{num}(:,:,t) = BW;
%         
%         if sum(idxOut) > length(idxOut)/2
%             break
%         end
%                 
%         figure(2); imshow( Is(:,:,t) ); title(num2str(t-1)); drawnow;
%         hold on;
%         plot(centroid(1), centroid(2), 'rx');
%         contour(BW, [.5 .5], 'r');
%         hold off;
%         drawnow;
%     end
% end

%%
if SAVE_VIDEO
    col = [ 0    0.4470    0.7410;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0.3010    0.7450    0.9330;...
            1.0000    0.0000    0.0000];

    outputVideo = VideoWriter([savefile '_soma.avi']);
    open(outputVideo);
    figure(4); set(gcf, 'Position', [0 0 1000 500]);
    for t = 0:num_timestep
        imshow( Is(:,:,t+1), []);
        set(gcf, 'Position', [0 0 1000 500]);
        hold on;
        for num = 1:length(flist)
            if nnz(soma{num}(:,:,t+1)) == 0, continue, end;
            contour(soma{num}(:,:,t+1), [.5 .5], 'Color', col(num,:), 'LineWidth', 2);
        end
        hold off;
        title(t);
        drawnow;
        writeVideo(outputVideo, getframe(gcf));
    end
    close(outputVideo);
end
disp('Done');
