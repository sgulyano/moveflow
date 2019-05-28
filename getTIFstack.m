addpath(genpath('toolbox'))

% filename = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/larva3_2_z%d_t%03d.tif';
% savefile = 'larva3_tif/larva3';
% num_timestep = 409;
% num_stacks = 8;
% imadj_param = [0.02,0.25];

filename = '~/Desktop/LipingsData/GFPRFPXYZTdata/larva4S14Zgreen/Larva 4_Series014_Crop001_t%03d_z%d_ch00.tif';
savefile = 'larva4s014_tif/larva4s014';
num_stacks = 2; 
num_timestep = 292;
imadj_param = [0.04 0.5];
        
for t = 0:num_timestep
    tifname = sprintf('%s_t%03d.tif', savefile, t);
    for z = 0:num_stacks
%         X = imread( sprintf(filename, z, t) );
        X = imread( sprintf(filename, t, z) );
        X = imadjust(X, imadj_param, [0 1]);
        if z == 0
            imwrite(X, tifname);
        else
            imwrite(X, tifname, 'WriteMode', 'append');
        end
    end
end
