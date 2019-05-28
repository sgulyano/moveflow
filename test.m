addpath(genpath('toolbox/utils'))
addpath(genpath('fix_slice'))

load('matlab.mat','directory', 'filename', 'aligndir', 'num_img', 'num_slice', 'switchTZ');
tic;
fix_slice_func(directory, filename, aligndir, num_img, num_slice, switchTZ);
toc;