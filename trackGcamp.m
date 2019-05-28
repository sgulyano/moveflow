close all
clear all
clc

addpath(genpath('toolbox'));
addpath(genpath('fix_slice'))
    
%% user parameters
DEBUG = true;

dataset = 3;
switchTZ = true;
soma_opt = struct('radii',[4 6], 'sensitiv',.95, 'edge_thr',.1);
stat = 'Mean';
switch dataset
    case 1
        % wild-type
        filepat = 'C:\Users\y_sar\Documents\MATLAB\DanData\GcamP images\WT\0322L1_100_P1\0322L1_100_P1_z%01d_t%03d.tif';
        savedir = 'WT0322L1_100_P1';
    case 2
        filepat = 'C:\Users\y_sar\Documents\MATLAB\DanData\GcamP images\WT\1110DLa_2_P1\1110DLa_2_P1_z%01d_t%02d.tif';
        savedir = 'WT1110DLa_2_P1';
    case 3
        % knockout
        filepat = 'C:\Users\y_sar\Documents\MATLAB\DanData\GcamP images\tmc1 mutant\20171213-L4-F_P1\20171213-L4-F_P1_z%02d_t%03d.tif';
        savedir = 'KO20171213-L4-F_P1';
end

col = lines(10);

%%
[di,ff,ext] = fileparts(filepat);
fl = dir(di);
valcell = arrayfun(@(x)sscanf(x.name, [ff ext]), fl, 'UniformOutput', false);
val = max([valcell{:}], [], 2);
if switchTZ
    n_sli = val(1);
    n_img = val(2);
else
    n_img = val(1);
    n_sli = val(2);
end

%% read middle slice for soma tracking
Vinfo = imfinfo(fullfile(di, [sprintf(ff,0,0) ext]));
V = zeros(Vinfo.Height, Vinfo.Width, n_img+1, 'uint8');
z = floor(n_sli/2);
for t = 0:n_img
    if switchTZ
        X = imread( fullfile(di, [sprintf(ff,z,t) ext]) );
    else
        X = imread( fullfile(di, [sprintf(ff,t,z) ext]) );
    end
    % normalize
    X = single(X);
    X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
    % adjust intensity
    n = 3;
    Idouble = im2double(X);
    avg = mean2(Idouble);
    sigma = std2(Idouble);
    X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
    V(:,:,t+1) = X;
end

%% let user select dendrite bounding box and circle soma
if exist(fullfile(savedir, [savedir '.mat']), 'file') ~= 2
    [ddad_soma, ddad_dend, ddae_soma, ddae_dend] = gui_select_neuron( V(:,:,1) );
    if exist(savedir, 'dir') ~= 7
        mkdir(savedir);
    end
    issave = true;
    if exist(fullfile(savedir, [savedir '.mat']), 'file') == 2
        save(fullfile(savedir, [savedir '.mat']), 'filepat', 'savedir', 'soma_opt', 'ddad_soma', 'ddad_dend', 'ddae_soma', 'ddae_dend', '-append');
    else
        save(fullfile(savedir, [savedir '.mat']), 'filepat', 'savedir', 'soma_opt', 'ddad_soma', 'ddad_dend', 'ddae_soma', 'ddae_dend');
    end
end

%% track all soma
load(fullfile(savedir, [savedir '.mat']), 'ddad_soma', 'ddad_dend', 'ddae_soma', 'ddae_dend');

listOfVariables = who('-file', fullfile(savedir, [savedir '.mat']));
if ismember('somas', listOfVariables)
    load(fullfile(savedir, [savedir '.mat']), 'somas');
else
    num_neuron = size(ddad_soma,1);

    centers = cell(1,n_img+1);
    radiis = cell(1,n_img+1);
    for t = 1:n_img+1
        I = imgaussian(double(V(:,:,t)),0.5);
        % for each frame
        [center, radii, metric] = imfindcircles(I, soma_opt.radii, ...
                'ObjectPolarity', 'bright', ...
                'Sensitivity', soma_opt.sensitiv, ...
                'EdgeThreshold',soma_opt.edge_thr, ...
                'Method','TwoStage');
        idx = metric > 0.2;
        centers{t} = center(idx,:);
        radiis{t} = radii(idx);
    end

    somas = soma_manual_multi( V, 1, [ddad_soma; ddae_soma], centers, radiis );
    save(fullfile(savedir, [savedir '.mat']), 'somas', '-append');

    ddad_somas = somas(1:num_neuron,:,:);
    ddae_somas = somas(num_neuron+1:end,:,:);

    % sanity check
    for t = 1:n_img+1
        figure(1); imshow(V(:,:,t));
        title(t-1);
        for i = 1:size(ddad_somas,1)
            viscircles(ddad_somas(i,1:2,t),ddad_somas(i,3,t), 'Color', col(i,:));
            text(ddad_somas(i,1,t),ddad_somas(i,2,t),'D','Color','r','FontSize',14);
            viscircles(ddae_somas(i,1:2,t),ddae_somas(i,3,t), 'Color', col(i,:));
            text(ddae_somas(i,1,t),ddae_somas(i,2,t),'E','Color','r','FontSize',14);
        end
        pause(0.1);
    end
end


%% track dendrites (as regions around soma)
listOfVariables = who('-file', fullfile(savedir, [savedir '.mat']));
if ismember('dds', listOfVariables)
    load(fullfile(savedir, [savedir '.mat']), 'dds');
else
    dds = dendrite_manual_multi( V, 1, [ddad_dend; ddae_dend], somas );
    save(fullfile(savedir, [savedir '.mat']), 'dds', '-append');

    num_neuron = size(ddad_dend,1);

    % sanity check
    ddad_dends = dds(1:num_neuron,:,:);
    ddae_dends = dds(num_neuron+1:end,:,:);
    for t = 1:n_img+1
        f = figure(2);
        set(f, 'Position', [100 100 1080 720]);
        imshow(V(:,:,t));
        title(t-1);
        hold on;
        for n = 1:num_neuron
            plot(ddad_dends(n,[1 3 3 1 1],t), ddad_dends(n,[2 2 4 4 2],t), 'g--', 'LineWidth', 2);
            text(mean(ddad_dends(n,[1 3],t)), mean(ddad_dends(n,[2 4],t)), 'D', 'FontSize', 14, 'Color', 'r');
            plot(ddae_dends(n,[1 3 3 1 1],t), ddae_dends(n,[2 2 4 4 2],t), 'g--', 'LineWidth', 2);
            text(mean(ddae_dends(n,[1 3],t)), mean(ddae_dends(n,[2 4],t)), 'E', 'FontSize', 14, 'Color', 'r');
        end
        hold off;
        pause(0.1);
    end
end

% % [HACK] swap ddaD and ddaE
% num_neuron = size(ddad_dend,1);
% tmp = ddad_dend; ddad_dend = ddae_dend; ddae_dend = tmp;
% tmp = ddad_soma; ddad_soma = ddae_soma; ddae_soma = tmp;
% dds = [dds(num_neuron+1:end,:,:); dds(1:num_neuron,:,:)];
% somas = [somas(num_neuron+1:end,:,:); somas(1:num_neuron,:,:)];
% save(fullfile(savedir, [savedir '.mat']), 'filepat', 'savedir', 'soma_opt', 'ddad_soma', 'ddad_dend', 'ddae_soma', 'ddae_dend', 'somas', 'dds');

%% calculate the intensity of soma and dendrites
[xx, yy] = meshgrid(1:size(V,2), 1:size(V,1));
n_neuron = size(somas,1);
dd_int = zeros(n_neuron, n_img+1);                  % dendrite intensity
sm_int = zeros(n_neuron, n_img+1);                  % soma intensity
dd_bws = zeros([size(V,1) size(V,2) n_img+1]);      % dendrite mask
sm_bws = zeros([size(V,1) size(V,2) n_img+1]);      % soma mask

nbins = 30;
inthist = zeros(n_neuron, n_img+1, nbins);
histedges = [0 linspace(4, 32, nbins-1) 255];

figure(3);
Itmp = zeros([size(V,1) size(V,2) 3]);
himg = imshow(Itmp);
hold on
hneuron = cell(1,n_neuron);
for n = 1:n_neuron
    Itmp(:,:,1) = col(n,1);
    Itmp(:,:,2) = col(n,2);
    Itmp(:,:,3) = col(n,3);
    hneuron{n} = imshow(Itmp);
end
hold off
for t = 1:n_img+1
    % for each frame
    Vt = zeros(Vinfo.Height, Vinfo.Width, n_sli+1, 'uint8');
    BW = zeros(Vinfo.Height, Vinfo.Width, n_sli+1, 'logical');
    for z = 0:n_sli
        % for each slice
        if switchTZ
            X = imread( fullfile(di, [sprintf(ff,z,t-1) ext]) );
        else
            X = imread( fullfile(di, [sprintf(ff,t-1,z) ext]) );
        end
        Vt(:,:,z+1) = X;
        
        % get mask
        I_thr = imgaussian(im2double(X), 1)>.1;
        I_erode = imerode(I_thr, strel('disk', 2));
        I_small = bwareaopen(I_erode, 20, 4);
        BW(:,:,z+1) = imdilate(I_small, strel('disk', 3));
    end
    It = max(Vt,[],3);
    BW2d = max(BW,[],3);
    
%     figure(5), imshowpair(It,BW2d);
%     set(5, 'Position', [100 100 1080 720]);
%     title(t-1);
%     pause(.1);

    
    dd_bw = zeros(size(It));
    sm_bw = zeros(size(It));
    for n = 1:n_neuron
        % for each soma
        sm_bw((xx - somas(n,1,t)).^2 + (yy - somas(n,2,t)).^2 <= (somas(n,3,t)+1)^2) = n;
        % for each dendrite
        mask = poly2mask(dds(n,[1 3 3 1 1],t), dds(n,[2 2 4 4 2],t), size(It,1), size(It,2));
        dd_bw(mask & ~BW2d) = n;
        % get intensity
        if strcmp(stat, 'Sum')
            sm_int(n,t) = sum(It(sm_bw == n));
            dd_int(n,t) = sum(It(dd_bw == n));
        elseif strcmp(stat, 'Mean')
            sm_int(n,t) = mean(It(sm_bw == n));
            dd_int(n,t) = mean(It(dd_bw == n));
        else
            error('Unknown Stat');
        end
        
        inthist(n,t,:) = histcounts(It(dd_bw == n), histedges);
        inthist(n,t,:) = inthist(n,t,:) / sum(inthist(n,t,:));
    end
    sm_bws(:,:,t) = sm_bw;
    dd_bws(:,:,t) = dd_bw;
    
    title(t-1);
    himg.CData = repmat(It, [1 1 3]);
    for n = 1:n_neuron
        hneuron{n}.AlphaData = double(dd_bws(:,:,t) == n) * 0.2;
    end
    pause(0.1);
end
%%
num_neuron = size(ddad_dend,1);
figure(30);
for ii = 1:num_neuron
    subplot(2,num_neuron,ii);
    imagesc(squeeze(inthist(ii,:,:)));
    xlabel('Bins');
    ylabel('Frame');
    title(['ddaD ' num2str(ii)]);
    subplot(2,num_neuron,ii+num_neuron);
    imagesc(squeeze(inthist(ii+num_neuron,:,:)));
    title(['ddaE ' num2str(ii)]);
    xlabel('Bins');
    ylabel('Frame');
end

%% visualize results
Vori = readImgSeq(filepat, n_img, n_sli);
%%
num_neuron = size(ddad_dend,1);

ddad_somas = somas(1:num_neuron,:,:);
ddae_somas = somas(num_neuron+1:end,:,:);
% sanity check
ddad_dends = dds(1:num_neuron,:,:);
ddae_dends = dds(num_neuron+1:end,:,:);

hdd = cell(num_neuron,2);
hgg = cell(num_neuron,1);
htt = cell(num_neuron,2);
hsm = cell(size(somas,1),size(dds,1));

outputVideo = VideoWriter(fullfile(savedir, [savedir '_' lower(stat) '.mp4']), 'MPEG-4');
open(outputVideo);

Itmp = zeros([size(V,1) size(V,2) 3]);

f = figure('Position', [100 0 1080 720], 'Name',savedir); clf; 
subplot(3,1,1); himg = imshow(Vori(:,:,1)); ttl = title('Time : 0');
hold on;
for n = 1:num_neuron
    Itmp(:,:,1) = col(n,1);
    Itmp(:,:,2) = col(n,2);
    Itmp(:,:,3) = col(n,3);
    hneuron{n} = imshow(Itmp);
    hneuron{n}.AlphaData = double(dd_bws(:,:,1) == n | dd_bws(:,:,1) == n+num_neuron) * 0.2;
end
for n = 1:num_neuron
    hdd{n,1} = plot(ddad_dends(n,[1 3 3 1 1],1), ddad_dends(n,[2 2 4 4 2],1), 'g--');
    hdd{n,2} = plot(ddae_dends(n,[1 3 3 1 1],1), ddae_dends(n,[2 2 4 4 2],1), 'c--');
    
    htt{n,1} = text(ddad_somas(n,1,1), ddad_somas(n,2,1), 'D', 'Color' ,'r', 'FontSize', 14);
    htt{n,2} = text(ddae_somas(n,1,1), ddae_somas(n,2,1), 'E', 'Color' ,'r', 'FontSize', 14);
    
    hgg{n} = viscircles([ddad_somas(n,1:2,1); ddae_somas(n,1:2,1)], [ddad_somas(n,3,1); ddae_somas(n,3,1)], 'EdgeColor', col(n,:));
end
hold off

hl = cell(1,4);
subplot(3,2,3); plot(0:n_img, sm_int(1:num_neuron,:)'); title(['ddaD Soma ' stat ' Intensity']);
hold on; hl{1} = plot([0 0], ylim, 'r', 'LineWidth', 2); hold off;
xlim([0 n_img]);
subplot(3,2,4); plot(0:n_img, sm_int(num_neuron+1:end,:)'); title(['ddaE Soma ' stat ' Intensity']);
hold on; hl{2} = plot([0 0], ylim, 'r', 'LineWidth', 2); hold off;
xlim([0 n_img]);

subplot(3,2,5); plot(0:n_img, dd_int(1:num_neuron,:)'); title(['ddaD Dendrite ' stat ' Intensity']);
hold on; hl{3} = plot([0 0], ylim, 'r', 'LineWidth', 2); hold off;
xlim([0 n_img]);
subplot(3,2,6); plot(0:n_img, dd_int(num_neuron+1:end,:)'); title(['ddaE Dendrite ' stat ' Intensity']);
hold on; hl{4} = plot([0 0], ylim, 'r', 'LineWidth', 2); hold off;
xlim([0 n_img]);

for t = 1:n_img+1
    ttl.String = ['Time : ' num2str(t-1)];

    subplot(3,1,1);
    himg.CData = Vori(:,:,t);
    for n = 1:num_neuron
        hdd{n,1}.XData = ddad_dends(n,[1 3 3 1 1],t);
        hdd{n,1}.YData = ddad_dends(n,[2 2 4 4 2],t);
        hdd{n,2}.XData = ddae_dends(n,[1 3 3 1 1],t);
        hdd{n,2}.YData = ddae_dends(n,[2 2 4 4 2],t);
        
        htt{n,1}.Position(1:2) = ddad_somas(n,1:2,t);
        htt{n,2}.Position(1:2) = ddae_somas(n,1:2,t);
        
        hneuron{n}.AlphaData = double(dd_bws(:,:,t) == n | dd_bws(:,:,t) == n+num_neuron) * 0.2;

        delete(hgg{n});
        hgg{n} = viscircles([ddad_somas(n,1:2,t); ddae_somas(n,1:2,t)], [ddad_somas(n,3,t); ddae_somas(n,3,t)], 'EdgeColor', col(n,:));
    end
    for i = 1:4
        hl{i}.XData = [t t];
    end
    pause(0.1);
    writeVideo(outputVideo, getframe(f));
end
close(outputVideo);

%% save to excel
T = table((0:n_img)');
T.Properties.VariableNames{1} = 'ImageNumber';

neu_row = cellstr([strcat('ddaD', cellfun(@num2str, num2cell(1:num_neuron)')); ...
        strcat('ddaE', cellfun(@num2str, num2cell(1:num_neuron)'))]);

soma_col = {'X', 'Y', 'R', 'Intensity'};
soma_tab = cat(2, somas, reshape(sm_int, [n_neuron 1 n_img+1]));
soma_tab = permute(soma_tab, [2 1 3]);
soma_tab = reshape(soma_tab, [n_neuron*4 n_img+1])';

Tsoma = array2table(soma_tab);
names = cellfun(@(x)strcat(x,'_',soma_col), neu_row, 'UniformOutput', false);
Tsoma.Properties.VariableNames = [names{:}];

dd_col = {'Left', 'Top', 'Right', 'Bottom', 'Intensity'};
dd_tab = cat(2, dds, reshape(dd_int, [n_neuron 1 n_img+1]));
dd_tab = permute(dd_tab, [2 1 3]);
dd_tab = reshape(dd_tab, [n_neuron*5 n_img+1])';

Tdend = array2table(dd_tab);
names = cellfun(@(x)strcat(x,'_',dd_col), neu_row, 'UniformOutput', false);
Tdend.Properties.VariableNames = [names{:}];

writetable([T Tsoma], fullfile(savedir, [savedir '_' lower(stat) '.xlsx']), 'Sheet', 1);
writetable([T Tdend], fullfile(savedir, [savedir '_' lower(stat) '.xlsx']), 'Sheet', 2);
