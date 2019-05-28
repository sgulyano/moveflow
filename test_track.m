close all; clc;

ROW = 1;
% filepath where vaa3d is installed
vaa3d = '~/Downloads/Vaa3d_V3.20_MacOSX10.9_64bit/Vaa3d.app/Contents/MacOS/vaa3d';
switch(1)
    case 1
        testmat = {'larva3_frames/larva3_track_fullflow.mat', ...               % list of traces to evaluate
                'larva3_frames/larva3_track.mat'};
        savefile = {'swc/fullflow/larva3_fullflow_n%d_t%03d.swc', ...
                'swc/our/larva3_our_n%d_t%03d.swc'};
        golddir = '~/Desktop/LipingsData/GoldStandard/ZstackL1_3_2grayscale/';  % directory containing the gold standard
        ofmat = 'larva3_frames/larva3_fullflow.mat';                            % Optical flow results
    case 2
        testmat = {'Larva4s014/Larva4s014_track_fullflow.mat', ... 
                'Larva4s014/Larva4s014_track.mat'};
        savefile = {'swc/fullflow/larva4S14Z_fullflow_n%d_t%03d.swc', ...
                'swc/our/larva4S14Z_our_n%d_t%03d.swc'};
        golddir = '~/Desktop/LipingsData/GoldStandard/larva4S14Z/';
        ofmat = 'Larva4s014/Larva4s014_fullflow.mat';                           % Optical flow results
end

%% compute SD, SSD, %SSD using vaa3d
addpath(genpath('toolbox'));

swcname = dir([golddir '*.swc']);

neuron_num = arrayfun(@(x)regexp(x.name,'.*_n(\d+)_t\d+.swc','tokens'), swcname);
neuron_num = cellfun(@(x)str2double(x{1}), neuron_num);

swc_time = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swcname);
swc_time = cellfun(@(x)str2double(x{1}), swc_time);

scores = cell(length(testmat),1);
for testnum = 1:length(testmat)
    disp(['Testing... ' testmat{testnum}]);
    load(testmat{testnum}, 'traces', 'initswctime');
    %score description (zero is good for all)
    %col1: entire-structure-average (from neuron 1 to 2)
    %col2: entire-structure-average (from neuron 2 to 1)
    %col3: average of bi-directional entire-structure-averages
    %col4: differen-structure-average
    %col5: percent of different-structure
    score = zeros(size(neuron_num,1),5);
    parfor i = 1:length(neuron_num)
        num = neuron_num(i);
    
        if swc_time(i) == initswctime(num), continue, end;
        swcgold = swcname(i).name;
        swctest = traces{swc_time(i)+1,num};

        % run neuron_distance in vaa3d
        goldfile = [golddir swcgold];

        testfile = sprintf(savefile{testnum}, num, swc_time(i));
        [directory, ~, ~] = fileparts(testfile);
        if ~exist(directory, 'dir'), mkdir(directory); end
        fileID = fopen(testfile, 'w');
        fprintf(fileID,'%d %d %.3f %.3f %.3f %.4f %d\n',swctest');
        fclose(fileID);

        outputfile = sprintf('output_n%d_t%03d.txt',num, swc_time(i));

        bashfile = sprintf('vaa3d_n%d_t%03d.sh',num, swc_time(i));
        fileID = fopen(bashfile, 'w');
        fprintf(fileID,'unset DYLD_FRAMEWORK_PATH DYLD_LIBRARY_PATH\n');
        fprintf(fileID,'%s -x neuron_distance -f neuron_distance -i %s %s -o %s\n', vaa3d, goldfile, testfile, outputfile);
        fclose(fileID);
        [status,cmdout] = system(['bash ' bashfile]);
        if status ~=0,
            disp(cmdout);
            error('Error: cannot compute neuron distance');
        end

        % read output file from vaa3d
        fileID = fopen(outputfile, 'r');
        fgetl(fileID);
        fgetl(fileID);
        for j = 1:5
            line = fgetl(fileID);
            token = regexp(line,'.* = (.*)','tokens');
            score(i,j) = str2double(token{1});
        end
    end
    scores{testnum} = score;
end
delete('output_n*_t*.txt')
delete('vaa3d_n*_t*.sh')
% delete(outputfile, testfile, outputfile);

%% sanity check
scores = cell(length(testmat),1);
for testnum = 1:length(testmat)
    disp(['Checking... ' testmat{testnum}]);
    for i = 1:length(neuron_num)
        num = neuron_num(i);
        if swc_time(i) == initswctime(num), continue, end;
        testfile = sprintf(savefile{testnum}, num, swc_time(i));
        data = readswc(testfile);
        figure(1); clf; EV_plot_img([], data); axis([0 512 0 256]);
        title(testmat{testnum});
        drawnow; pause(0.1);
    end
end

%% plot results framewise
figure; set(gcf,'color','w','position',[0 0 800 200]);
for nnum = unique(neuron_num)'
    fig = mod(nnum-1,ROW);
    if nnum > ROW, figure; set(gcf,'color','w','position',[0 0 800 200]); end;
    framestr = arrayfun(@num2str,swc_time(neuron_num==nnum),'UniformOutput',false);
    swc_ti = swc_time(neuron_num==nnum);
    idx = diff(swc_ti,1)>1;
    idx = [1; idx]; idx(end) = 1;
    for k = 1:length(idx)
        if idx(k) == 0, framestr{k} = ''; end
    end
    
    val = cellfun(@(x)(x(neuron_num==nnum,3)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+1); plot([val{:}], 'LineWidth', 2)
    ylabel('SD');
    xlabel('frames');
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr);
    ylim([0 15])
    h_legend = legend({'FullFlow','Ours'}, 'location', 'northwest');
    set(h_legend,'FontSize',14);

    val = cellfun(@(x)(x(neuron_num==nnum,4)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+2); plot([val{:}], 'LineWidth', 2)
    ylabel('SSD');
    xlabel('frames');
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr);
    title(['Neuron ' num2str(nnum)]);
    ylim([0 15])
    
    val = cellfun(@(x)(x(neuron_num==nnum,5)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+3); plot([val{:}], 'LineWidth', 2)
    ylabel('%SSD');
    xlabel('frames');
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr);
    ylim([0 1])
end

%% plot box plot
figure; set(gcf,'color','w','position',[0 0 250 250]);
for nnum = unique(neuron_num)'
    fig = mod(nnum-1,ROW);
    if nnum > ROW, figure; set(gcf,'color','w','position',[0 0 250 250]); end;
    framestr = arrayfun(@num2str,swc_time(neuron_num==nnum),'UniformOutput',false);
    swc_ti = swc_time(neuron_num==nnum);
    idx = diff(swc_ti,1)>1;
    idx = [1; idx]; idx(end) = 1;
    for k = 1:length(idx)
        if idx(k) == 0, framestr{k} = ''; end
    end
    
    val = cellfun(@(x)(x(neuron_num==nnum,3)), scores, 'UniformOutput', false');
    sd = [val{:}];
    boxplot(sd);
    title(['Neuron ' num2str(nnum)]);
    ylabel('SD (pixels)');
    set(gca,'xtick',1:2);
    set(gca,'xticklabel',{'FullFlow','Ours'});
end

%% plot trace in each frame
load(ofmat, 'gfpadj');
A = cell(length(testmat),1);
for i = 1:length(testmat);
    A{i} = load(testmat{i}, 'traces');
end

for i = 3:length(neuron_num)
    num = neuron_num(i);
    t = swc_time(i);
    
    swcgold = swcname(i).name;
    goldtrace = read_swc_file([golddir swcgold]);
    
    figure(99), set(gcf, 'color', 'w', 'position', [100 100 800 400]);
    I = im2double(gfpadj(:,:,t+1));
    MAXDIST = 10;
    ii = round(goldtrace(:,4)); jj = round(goldtrace(:,3));
    xmin = max(min(jj)-MAXDIST*2, 1); ymin = max(min(ii)-MAXDIST*2, 1);
    xmax = min(max(jj)+MAXDIST, size(I,2)); ymax = min(max(ii)+MAXDIST, size(I,1));
    I = I(ymin:ymax, xmin:xmax);
    
    subplot(1,3,1); imshow(I); title(t);
    hold on; EV_plot_img( [], bsxfun(@minus, goldtrace, [0 0 xmin-1 ymin-1 0 0 0]), 'r-o' );
    subplot(1,3,2); imshow(I); title('FullFlow');
    hold on; EV_plot_img( [], bsxfun(@minus, A{1}.traces{t+1,num}, [0 0 xmin-1 ymin-1 0 0 0]), 'r-o' );
    subplot(1,3,3); imshow(I); title('Ours');
    hold on; EV_plot_img( [], bsxfun(@minus, A{2}.traces{t+1,num}, [0 0 xmin-1 ymin-1 0 0 0]), 'r-o' );
    pause;
end


%%

%% plot results framewise
figure; set(gcf,'color','w','position',[0 0 800 200]);
for nnum = 3%unique(neuron_num)'
    
    framestr = arrayfun(@num2str,swc_time(neuron_num==nnum),'UniformOutput',false);
    swc_ti = swc_time(neuron_num==nnum);
    idx = diff(swc_ti,1)>1;
    idx = [1; idx]; idx(end) = 1;
    for k = 1:length(idx)
        if idx(k) == 0, framestr{k} = ''; end
    end
    
    subplot(1,3,1);
    val = cellfun(@(x)(x(neuron_num==nnum,3)), scores, 'UniformOutput', false');
    sd = [val{:}];
    bh=boxplot(sd);
    set(bh,'linewidth',2);
    
    ylabel('SD (pixels)');
    set(gca,'xtick',1:2);
    set(gca,'xticklabel',{'FullFlow','Ours'});
    
    
    framestr = arrayfun(@num2str,swc_time(neuron_num==nnum),'UniformOutput',false);
    swc_ti = swc_time(neuron_num==nnum);
    idx = diff(swc_ti,1)>1;
    idx = [1; idx]; idx(end) = 1;
    for k = 1:length(idx)
        if idx(k) == 0, framestr{k} = ''; end
    end

    val = cellfun(@(x)(x(neuron_num==nnum,4)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+2); plot([val{:}], 'LineWidth', 3)
    ylabel('SSD');
    xlabel('frames');
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr);
    title(['Neuron ' num2str(nnum)]);
    h_legend = legend({'FullFlow','Ours'}, 'location', 'northwest');
    set(h_legend,'FontSize',14);
    ylim([0 15])
    
    val = cellfun(@(x)(x(neuron_num==nnum,5)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+3); plot([val{:}], 'LineWidth', 3)
    ylabel('%SSD');
    xlabel('frames');
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr);
    ylim([0 1])
end