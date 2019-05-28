close all; clear all; clc;

% filepath where vaa3d is installed
vaa3d = '~/Downloads/Vaa3d_V3.20_MacOSX10.9_64bit/Vaa3d.app/Contents/MacOS/vaa3d';
switch(1)
    case 1
        golddir = '~/Desktop/LipingsData/GoldStandard/ZstackL1_3_2grayscale/';  % directory containing the gold standard
        gfpfile = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/larva3_2_z%d_t%03d.tif';     % Ca images
        switchTZ = true;                                % switch position of T and Z in filename
        normalized = true;                              % normalized image slice
        use_normxcorr = false;                          % use Normalized Cross-Correlation to track soma
        adj_slicewise = false;                          % apply intensity adjustment slice-wise
        num_slice = 8;                                  % number of image slices per stack
        num_img = 409;                                  % number of frames (start from 0 to num_img)
        gfpadjrange = [0.02 0.1];                      % intensity adjusting threshold for Ca images
        use_medfilt = false;                            % use median filter in soma detection               
        savefile = 'swc/neutube/larva3_neutube_n%d_t%03d.swc';
        savefile_one = 'swc/neutubeone/larva3_neutubeone_n%d_t%03d.swc';
    case 2
        golddir = '~/Desktop/LipingsData/GoldStandard/larva4S14Z/';
        gfpfile = '~/Desktop/LipingsData/GFPRFPXYZTdata/larva4S14Zgreen/Larva 4_Series014_Crop001_t%03d_z%d_ch00.tif';
        switchTZ = false;
        normalized = false;                             % normalized image slice
        use_normxcorr = false;                          % use Normalized Cross-Correlation to track soma
        adj_slicewise = false;                          % apply intensity adjustment slice-wise
        num_slice = 2;      % number of image slices per stack
        num_img = 292;      % number of frames (start from 0 to num_img)
        gfpadjrange = [0.04 0.2]; % for FullFlow use [0.04 0.3]
        use_medfilt = '3D';
        savefile = 'swc/neutube/larva4S14Z_neutube_n%d_t%03d.swc';
        savefile_one = 'swc/neutubeone/larva4S14Z_neutubeone_n%d_t%03d.swc';
end
        
%% compute SD, SSD, %SSD using vaa3d
addpath(genpath('toolbox'));

gfpinfo = imfinfo(sprintf(gfpfile,0,0));
swcname = dir([golddir '*.swc']);

neuron_num = arrayfun(@(x)regexp(x.name,'.*_n(\d+)_t\d+.swc','tokens'), swcname);
neuron_num = cellfun(@(x)str2double(x{1}), neuron_num);

swc_time = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swcname);
swc_time = cellfun(@(x)str2double(x{1}), swc_time);

%% trace using neuTube
%score description (zero is good for all)
%col1: entire-structure-average (from neuron 1 to 2)
%col2: entire-structure-average (from neuron 2 to 1)
%col3: average of bi-directional entire-structure-averages
%col4: differen-structure-average
%col5: percent of different-structure
score = zeros(size(neuron_num,1),5);
for i = 1:length(neuron_num)
    %% trace each frame independently using neuTube
    num = neuron_num(i);
    t = swc_time(i);
    
    fprintf('Time : %d\n', t);
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

    stackfile = 'stack.tif';
    for k = 1:size(V,3)
        V(:,:,k) = imadjust(V(:,:,k), gfpadjrange, [0 1]);
        if k == 1
            imwrite(flipud(V(:,:,1)), stackfile);
        else
            imwrite(flipud(V(:,:,k)), stackfile, 'writemode', 'append');
        end
    end

    bashfile = 'vaa3d_neutube.sh';
    fileID = fopen(bashfile, 'w');
    fprintf(fileID,'unset DYLD_FRAMEWORK_PATH DYLD_LIBRARY_PATH\n');
    fprintf(fileID,'%s -x neuTube -m tracing -i %s\n', vaa3d, stackfile);
    fclose(fileID);
    [status,cmdout] = system(['bash ' bashfile]);
    if status ~=0,
        disp(cmdout);
        error('Error: cannot compute neuron distance');
    end
    
    %% remove irrelevant neuron traces
    % get bounding box
    goldfile = [golddir swcname(i).name];
    goldswc = readswc(goldfile);
    box = [min(goldswc(:,3:4)) max(goldswc(:,3:4))];
    
    testswc = readswc( sprintf('%s_neutube.swc', stackfile) );
    figure(1); subplot(2,1,1); EV_plot_img( max(V,[],3), testswc );
    hold on; plot(box([1 3 3 1 1]), box([2 2 4 4 2]), 'g--'); hold off;
    title( sprintf('Time : %d', t) );
    
    ed = testswc(:, [1 7]);
    ed = ed(all(ed > 0,2),:);
    G = sparse(ed(:,1), ed(:,2), 1, size(testswc,1), size(testswc,1));
    G = G + G';
    [ci, sizes] = components(G);
    
    % keep only neurons inside the bounding box
    st_idx = 0;
    data = zeros(0,7);
    for k = 1:length(sizes)
        nc = sum(ci==k);
        idx = all(testswc(ci==k,3:4) >= ones(nc,1)*box(1:2) & testswc(ci==k,3:4) <= ones(nc,1)*box(3:4),2);
        if sum(idx) > 0            
            % convert to swc format
            datak = graph2swc( G(ci==k,ci==k), testswc(ci==k,3:6) );
            datak(:,1) = datak(:,1) + st_idx;
            datak(2:end,7) = datak(2:end,7) + st_idx;
            data = [data; datak];
            st_idx = st_idx + nc;
        end
    end
    
    figure(1); subplot(2,1,2); EV_plot_img( max(V,[],3), data );
    hold on; plot(box([1 3 3 1 1]), box([2 2 4 4 2]), 'g--'); hold off;
    drawnow;
    % write to file
    testfile = sprintf(savefile,num, t);
    [directory, ~, ~] = fileparts(testfile);
    if ~exist(directory, 'dir'), mkdir(directory); end
    fileID = fopen(testfile, 'w');
    fprintf(fileID,'%d %d %.3f %.3f %.3f %.4f %d\n',data');
    fclose(fileID);
end

delete 'stack.tif' 'stack.tif_neutube.swc' 'vaa3d_neutube.sh'

%% keep only one tree
for i = 1:length(neuron_num)
    num = neuron_num(i);
    t = swc_time(i);
    
    testswc = readswc( sprintf(savefile, num, t) );
    
    if isempty(testswc)
        copyfile(sprintf(savefile, num, t), sprintf(savefile_one, num, t));
        continue
    end
    
    figure(1); clf; subplot(2,1,1); EV_plot_img([], testswc );
    axis([0 512 0 256]);
    title(t);
    drawnow;
    
    ed = testswc(:, [1 7]);
    ed = ed(all(ed > 0,2),:);
    G = sparse(ed(:,1), ed(:,2), 1, size(testswc,1), size(testswc,1));
    G = G + G';
    [ci, sizes] = components(G);
    [~, pos] = max(sizes);
    % keep only the largest tree
    nc = sum(ci==pos);
    data = graph2swc( G(ci==pos,ci==pos), testswc(ci==pos,3:6) );

    
    figure(1); subplot(2,1,2); EV_plot_img([], data );
    axis([0 512 0 256]);
    drawnow;
    
    % write to file
    testfile = sprintf(savefile_one, num, t);
    [directory, ~, ~] = fileparts(testfile);
    if ~exist(directory, 'dir'), mkdir(directory); end
    fileID = fopen(testfile, 'w');
    fprintf(fileID,'%d %d %.3f %.3f %.3f %.4f %d\n',data');
    fclose(fileID);
    pause(.1);
end