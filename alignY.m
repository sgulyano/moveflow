close all; clear all; clc;

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
    case 6
        gfpfile = '~/Desktop/LipingsData/20170616L3/20170616L3_z%02d_t%03d.tif';     % Ca images
        rfpfile = '';
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 10;                                 % number of image slices per stack
        end_slice = 10;                                 % image slice to stop reading
        num_img = 711;                                  % number of frames (start from 0 to num_img)
        savefile = '20170616L3/20170616L3';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 72;                                % frame picked for manual tracing
        motion_end_frames = 111;                        % frame that motion end correspond to pickframes
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        soma_opt = struct();                            % options for soma detector
        soma_opt.radii = [15 30];
        soma_opt.sensitiv = 0.6;
        soma_opt.edge_thr = 0.3;
        adj_soma = [29, 86];
    case 7
        gfpfile = '~/Desktop/LipingsData/20170616L7/20170616L7_z%d_t%03d.tif';     % Ca images
        gfpfile_align = '~/Desktop/LipingsData/20170616L7_align/20170616L7_align_z%d_t%03d.tif';     % Ca images
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
        pickframes = 72;                                % frame picked for manual tracing
        motion_end_frames = 111;                        % frame that motion end correspond to pickframes
        pickslice = 5;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        soma_opt = struct();                            % options for soma detector
        soma_opt.radii = [15 30];
        soma_opt.sensitiv = 0.6;
        soma_opt.edge_thr = 0.3;
        adj_soma = [29, 86];
end

%% read img data
disp('Reading GFP and RFP images');
Vinfo = imfinfo(sprintf(gfpfile,0,0));

gfp = zeros(Vinfo.Height, Vinfo.Width, num_img+1, 'uint8');
gfp_align = zeros(Vinfo.Height, Vinfo.Width, num_img+1, 'uint8');
for t = 0:num_img
    V = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, 'uint8');
    Vori = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, 'uint8');
    for z = 0:num_slice
        if switchTZ
            Xori = imread( sprintf(gfpfile,z,t) );
        else
            Xori = imread( sprintf(gfpfile,t,z) );
        end
        X = Xori;
        if normalized
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        end
        
        if adj_slicewise
            n = 3;
            Idouble = im2double(X);
            avg = mean2(Idouble);
            sigma = std2(Idouble);
            X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
        end
        Vori(:,:,z+1) = Xori;
        V(:,:,z+1) = X;
    end
    gfp(:,:,t+1) = max(V,[],3);
    
    bs = padarray(V(:,:,pickslice), [10, 0], 'both');
    Valign = V;
    for z = 0:num_slice
        C = normxcorr2(V(:,:,z+1), bs);
        [ypeak, xpeak] = find(C==max(C(:)));
        yoffSet = ypeak-size(V,1);
        xoffSet = xpeak-size(V,2);

        Valign(:,:,z+1) = imtranslate(V(:,:,z+1),[0 yoffSet-10]);
        
        Xori = imtranslate(Vori(:,:,z+1),[0 yoffSet-10]);
        
        if switchTZ
            imwrite(Xori, sprintf(gfpfile_align,z,t));
        else
            imwrite(Xori, sprintf(gfpfile_align,t,z));
        end
    end
    
%     gfp_align(:,:,t+1) = max(Valign,[],3);
    
    fprintf('.');
    if mod(t+1,100) == 0, fprintf('\n'); end;
end
fprintf('\n');
figure(1), imshow3D(gfp)
% figure(2), imshow3D(gfp_align)

% %%
% bs = padarray(V(:,:,pickslice), [10, 0], 'both');
% Valign = V;
% for z = 0:num_slice
%     C = normxcorr2(V(:,:,z+1), bs);
%     [ypeak, xpeak] = find(C==max(C(:)));
%     yoffSet = ypeak-size(V,1);
%     xoffSet = xpeak-size(V,2);
%     
%     Valign(:,:,z+1) = imtranslate(V(:,:,z+1),[0 yoffSet-10]);
% end
% figure(2); subplot(2,1,1); imshow(max(V,[],3));
% subplot(2,1,2); imshow(max(Valign,[],3))


%%
outputVideo = VideoWriter([savefile '_Yalign.avi']);
open(outputVideo);
figure(2); clf;

disp('Read Aligned Data');
for t = 0:num_img
    V = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, 'uint8');
    Vori = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, 'uint8');
    for z = 0:num_slice
        if switchTZ
            Xori = imread( sprintf(gfpfile_align,z,t) );
        else
            Xori = imread( sprintf(gfpfile_align,t,z) );
        end
        X = Xori;
        if normalized
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        end
        
        if adj_slicewise
            n = 3;
            Idouble = im2double(X);
            avg = mean2(Idouble);
            sigma = std2(Idouble);
            X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
        end
        Vori(:,:,z+1) = Xori;
        V(:,:,z+1) = X;
    end
    gfp_align(:,:,t+1) = max(V,[],3);
    
    figure(2); imshow(max(V,[],3));
    title(t)
    set(gcf, 'Position', [0 500 800 500]);
    drawnow;
    pause(0.01);
    writeVideo(outputVideo, getframe(gcf));
    
    fprintf('.');
    if mod(t+1,100) == 0, fprintf('\n'); end;
end
figure(3); imshow3D(gfp_align);
close(outputVideo);
