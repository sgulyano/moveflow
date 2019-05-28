gfpfile = '~/Desktop/LipingsData/1122larva5_3_Subset for Sarun/1122larva5_3_Subset_z%d_t%03d.tif';     % Ca images
rfpfile = '';
switchTZ = true;                                % switch position of T and Z in filename
normalized = true;                              % normalized image slice
start_slice = 2;                                % image slice to start reading
num_slice = 5;                                  % number of image slices per stack
num_img = 256;                                  % number of frames (start from 0 to num_img)
savefile = '1122larva5_3/1122larva5_3';
gfpadjrange = [0.02,0.15];                      % intensity adjusting threshold for Ca images
pickslice = 4;                                  % image slice for tracking soma
use_medfilt = false;                            % use median filter in soma detection
soma_opt = struct();                            % options for soma detector

%% read img data
disp('Reading GFP and RFP images');
gfpinfo = imfinfo(sprintf(gfpfile,0,0));
gfp = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
for t = 161%0:num_img
    V = zeros(gfpinfo.Height, gfpinfo.Width, num_slice+1, 'uint8');
    for z = start_slice:num_slice
        if switchTZ
            X = imread( sprintf(gfpfile,z,t) );
        else
            X = imread( sprintf(gfpfile,t,z) );
        end
        if normalized
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        end
        V(:,:,z+1) = X;
%         if z == pickslice
%             gfp1sl(:,:,t+1) = X;
%         end
    end
%     gfp(:,:,t+1) = max(V,[],3);
    
%     if ~isempty(rfpfile)
%         V = zeros(rfpinfo.Height, rfpinfo.Width, num_slice+1);
%         for z = 0:num_slice
%             V(:,:,z+1) = imread( sprintf(rfpfile,t,z) );
%         end
%         rfp(:,:,t+1) = max(V,[],3);
%     end
%     fprintf('.');
%     if mod(t+1,100) == 0, fprintf('\n'); end;
end
%%
Vadj = V(:,:,start_slice+1:num_slice+1);
for z = 1:size(Vadj,3)
    n = 2;
    Idouble = im2double(Vadj(:,:,z));
    avg = mean2(Idouble);
    sigma = std2(Idouble);
    Vadj(:,:,z) = imadjust(Vadj(:,:,z),[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
end
figure, imshow3D(Vadj);