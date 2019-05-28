close all; clc;

addpath('../toolbox/utils');
addpath(genpath('nonrigid_version23'));

%% user parameters
dataset = 8;
switch dataset
    case 8
        gfpfile = '~/Desktop/LipingsData/20171228L2/20171228L2_z%02d_t%03d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 635;                                  % number of frames (start from 0 to num_img)
        savefile = '20171228L2/20171228L2';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = [1 397 462 511];                   % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = false;                              % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 9
        gfpfile = '~/Desktop/LipingsData/20171228L1_B/20171228L1_B_z%02d_t%03d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 12;                                 % number of image slices per stack
        end_slice = 12;                                 % image slice to stop reading
        num_img = 201;                                  % number of frames (start from 0 to num_img)
        savefile = '20171228L1_B/20171228L1_B';
        gfpadjrange = [0.15, 1];                        % intensity adjusting threshold for Ca images
        pickframes = 1;                                 % frame picked for manual tracing
        pickslice = 6;                                  % image slice for tracking soma
        use_medfilt = '2D';                             % use median filter in soma detection
        ddaD_Left = false;                              % Is ddaD on the left?
        soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
    case 10
        gfpfile = '~/Desktop/LipingsData/20171228L5/20171228L5_z%02d_t%03d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 328;                                  % number of frames (start from 0 to num_img)
        savefile = '20171228L5/20171228L5';
        pickslice = 6;                                  % image slice for tracking soma
    case 11
        gfpfile = '~/Desktop/LipingsData/20171229L2/20171229L2_z%02d_t%02d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 92;                                   % number of frames (start from 0 to num_img)
        savefile = '20171229L2/20171229L2';
        pickslice = 6;                                  % image slice for tracking soma
    case 12
        gfpfile = '~/Desktop/LipingsData/20171229L8-B/20171229L8_B_good_z%02d_t%03d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 13;                                 % number of image slices per stack
        end_slice = 13;                                 % image slice to stop reading
        num_img = 205;                                  % number of frames (start from 0 to num_img)
        savefile = '20171229L8_B/20171229L8_B';
        pickslice = 6;                                  % image slice for tracking soma
    case 13
        gfpfile = '~/Desktop/LipingsData/20171229L12_B/20171229L12_B_good_z%02d_t%02d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        num_slice = 13;                                 % number of image slices per stack
        num_img = 71;                                  % number of frames (start from 0 to num_img)
        savefile = '20171229L12_B/20171229L12_B';
        pickslice = 6;                                  % image slice for tracking soma
    case 14
        gfpfile = '~/Desktop/LipingsData/20171229L13/20171229L13_F_z%02d_t%02d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        adj_slicewise = true;                           % apply intensity adjustment slice-wise
        use_normxcorr = true;                           % use Normalized Cross-Correlation to track soma
        start_slice = 0;                                % image slice to start reading
        num_slice = 11;                                 % number of image slices per stack
        end_slice = 11;                                 % image slice to stop reading
        num_img = 106;                                  % number of frames (start from 0 to num_img)
        savefile = '20171229L13/20171229L13';
        pickslice = 5;                                  % image slice for tracking soma
end

%% read img data
disp('Reading GFP and RFP images');
gfpinfo = imfinfo(sprintf(gfpfile,0,0));
gfp = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
gfp2 = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,num_img) '\n\n']);

% figure(1); clf; himg = imshow(gfp(:,:,1));

parfor t = 0:num_img
    V = zeros(gfpinfo.Height, gfpinfo.Width, num_slice+1, 'uint8');
    for z = 0:num_slice
        if switchTZ
            X = imread( sprintf(gfpfile,z,t) );
        else
            X = imread( sprintf(gfpfile,t,z) );
        end
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
        V(:,:,z+1) = X;
    end
    gfp2(:,:,t+1) = max(V,[],3);
    
    %% align slices
    Valign = align_stack(V, pickslice+1);
    gfp(:,:,t+1) = uint8(max(Valign,[],3)*255);

    
%     set(himg, 'CData', gfp(:,:,t+1));
%     title(t); drawnow;
    fprintf('\b|\n');
end
fprintf('\n');

figure(2); clf; imshow3D(gfp); title(savefile);
figure(3); clf; imshow3D(gfp2); title(savefile);
save([fullfile('..',savefile) '_align_gfp.mat'], 'gfp');