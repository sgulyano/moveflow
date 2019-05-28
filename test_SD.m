close all; clc;

SAVEXLSX = false;

% filepath where vaa3d is installed
vaa3d = '~/Downloads/Vaa3d_V3.20_MacOSX10.9_64bit/Vaa3d.app/Contents/MacOS/vaa3d';

switch(1)
    case 1
        tracefile = 'swc/%s/larva3_%s_n%d_t%03d.swc';
        golddir = '~/Desktop/LipingsData/GoldStandard/ZstackL1_3_2grayscale/';  % directory containing the gold standard
        initswctime = [198, 201, 229, 252, 340, 342];
        xlsname = 'larva3_SSDdata.xlsx';
        ofmat = 'larva3_frames/larva3_fullflow.mat';                            % Optical flow results
    case 2
        tracefile = 'swc/%s/larva4S14Z_%s_n%d_t%03d.swc';
        golddir = '~/Desktop/LipingsData/GoldStandard/larva4S14Z/';
        initswctime = 1;
        xlsname = 'larva4S14Z_SSDdata.xlsx';
        ofmat = 'Larva4s014/Larva4s014_fullflow.mat';                           % Optical flow results
end
methods = {'neutube', 'fullflow', 'our'};
%% compute SD, SSD, %SSD using vaa3d
addpath(genpath('toolbox'));

swcname = dir([golddir '*.swc']);

neuron_num = arrayfun(@(x)regexp(x.name,'.*_n(\d+)_t\d+.swc','tokens'), swcname);
neuron_num = cellfun(@(x)str2double(x{1}), neuron_num);

swc_time = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swcname);
swc_time = cellfun(@(x)str2double(x{1}), swc_time);

scores = zeros(length(neuron_num), length(methods), 5);
parfor i = 1:length(neuron_num)
    num = neuron_num(i);
    t = swc_time(i);
    if any(t == initswctime), continue, end;
    goldfile = [golddir swcname(i).name];
    score = zeros(length(methods), 5);
    for j = 1:length(methods)
        %score description (zero is good for all)
        %col1: entire-structure-average (from neuron 1 to 2)
        %col2: entire-structure-average (from neuron 2 to 1)
        %col3: average of bi-directional entire-structure-averages
        %col4: differen-structure-average
        %col5: percent of different-structure
        testfile = sprintf(tracefile, methods{j}, methods{j}, num, t);
        
        testswc = readswc(testfile);
        if isempty(testswc), testfile = 'swc/swcempty.swc'; end
        
        bashfile = sprintf('vaa3d_n%d.sh', i);
        outputfile = sprintf('output_n%d.txt', i);
        
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
        for k = 1:5
            line = fgetl(fileID);
            token = regexp(line,'.* = (.*)','tokens');
            val = str2double(token{1});
            if val < 0
                keyboard;
            end
            score(j,k) = str2double(token{1});
        end
    end
    scores(i,:,:) = score;
end
delete('output_n*.txt')
delete('vaa3d_n*.sh')
disp('Done');
%%
if SAVEXLSX
    sheetname = {'structavg12', 'structavg21', 'SD', 'SSD', 'percentSSD'};
    for kk = 1:5
        T = array2table([swc_time, neuron_num, scores(:,:,kk)]);
        T.Properties.VariableNames = {'Frame','NeuronNumber',methods{:}};
        writetable(T, xlsname,'Sheet',sheetname{kk})
    end
end

%% sanity check
for testnum = 1:length(methods)
    disp(['Checking... ' methods{testnum}]);
    for i = 1:length(neuron_num)
        num = neuron_num(i);
        t = swc_time(i);
        if swc_time(i) == initswctime(num), continue, end;
        testfile = sprintf(tracefile, methods{testnum}, methods{testnum}, num, t);
        data = readswc(testfile);
        figure(1); clf; EV_plot_img([], data); axis([0 512 0 256]);
        title(methods{testnum});
        drawnow; pause(0.1);
    end
end

%% plot results framewise
ROW = 1;
gtidx = all(scores(:,:,3)==0,2);
figure; set(gcf,'color','w','position',[0 0 1400 250]);
for nnum = 1:4%unique(neuron_num)'
    fig = mod(nnum-1,ROW);
    if nnum > ROW, figure; set(gcf,'color','w','position',[0 0 1400 250]); end;
    numi = neuron_num==nnum & ~gtidx;
    framestr = arrayfun(@num2str,swc_time(numi),'UniformOutput',false);
    swc_ti = swc_time(numi);
    idx = diff(swc_ti,1)>3;
    idx = [1; idx]; idx(end) = 1;
    for k = 1:length(idx)
        if idx(k) == 0, framestr{k} = ''; end
    end

    val = scores(numi,:,3); %cellfun(@(x)(x(neuron_num==nnum,3)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+1); plot(val, 'LineWidth', 3)
    ylabel('SD', 'FontSize',14);
    xlabel('frames', 'FontSize',14);
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr, 'FontSize',14);
    axis tight
    ylim([0 30])
    h_legend = legend({'neuTube', 'FullFlow','Ours'}, 'location', 'northwest');
    set(h_legend,'FontSize',14);

    val = scores(numi,:,4); %cellfun(@(x)(x(neuron_num==nnum,4)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+2); plot(val, 'LineWidth', 3)
    ylabel('SSD', 'FontSize',14);
    xlabel('frames', 'FontSize',14);
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr, 'FontSize',14);
    title(['Neuron ' num2str(nnum)]);
    axis tight
    ylim([0 30])
    h_legend = legend({'neuTube', 'FullFlow','Ours'}, 'location', 'northwest');
    set(h_legend,'FontSize',14);
    
    val = scores(numi,:,5); %cellfun(@(x)(x(neuron_num==nnum,5)), scores, 'UniformOutput', false');
    subplot(ROW,3,fig*3+3); plot(val, 'LineWidth', 3)
    ylabel('%SSD', 'FontSize',14);
    xlabel('frames', 'FontSize',14);
    set(gca,'xtick',1:length(framestr));
    set(gca,'xticklabel',framestr, 'FontSize',14);
    axis tight
    ylim([0 1])
end

%% plot box plot
figure; set(gcf,'color','w','position',[0 0 400 250]);
for nnum = 1:4%unique(neuron_num)'
    fig = mod(nnum-1,ROW);
    if nnum > ROW, figure; set(gcf,'color','w','position',[0 0 400 250]); end;
    numi = neuron_num==nnum & ~gtidx;
    framestr = arrayfun(@num2str,swc_time(numi),'UniformOutput',false);
    swc_ti = swc_time(numi);
    idx = diff(swc_ti,1)>3;
    idx = [1; idx]; idx(end) = 1;
    for k = 1:length(idx)
        if idx(k) == 0, framestr{k} = ''; end
    end
    val = scores(numi,:,3);
    sd = val;
    h = boxplot(sd, 'Orientation', 'horizontal');
    set(h,{'linew'},{2})
    title(['Neuron ' num2str(nnum)]);
    xlabel('SD (pixels)', 'FontSize',14);
    set(gca,'ytick',1:3);
    set(gca,'yticklabel',{'neuTube', 'FullFlow','Ours'}, 'FontSize',14);
    xlim([0 30])
end

%% plot trace in each frame
load(ofmat, 'gfpadj');
methods{1} = 'neutubeone';
for i = 1:length(neuron_num)
    num = neuron_num(i);
    t = swc_time(i);
    if swc_time(i) == initswctime(num), continue, end;
    swcgold = swcname(i).name;
    goldtrace = read_swc_file([golddir swcgold]);
    
    figure(99), set(gcf, 'color', 'w', 'position', [100 100 1200 400]);
    I = im2double(gfpadj(:,:,t+1));
    MAXDIST = 10;
    ii = round(goldtrace(:,4)); jj = round(goldtrace(:,3));
    xmin = max(min(jj)-MAXDIST*2, 1); ymin = max(min(ii)-MAXDIST*2, 1);
    xmax = min(max(jj)+MAXDIST, size(I,2)); ymax = min(max(ii)+MAXDIST, size(I,1));
    I = I(ymin:ymax, xmin:xmax);
    
    subplot(1,4,1); imshow(I); title(t);
    hold on;
    EV_plot_img( [], bsxfun(@minus, goldtrace, [0 0 xmin-1 ymin-1 0 0 0]), 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', 2,'markers',8);
    hold off;    
    for j = 1:3
        swcfile = sprintf(tracefile, methods{j}, methods{j}, num, t);
        traces = readswc(swcfile);
        if isempty(traces), continue, end
        subplot(1,4,j+1); imshow(I); title(methods{j});
        hold on; 
        EV_plot_img( [], bsxfun(@minus, traces, [0 0 xmin-1 ymin-1 0 0 0]), 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', 2,'markers',8);
        hold off;
    end
%     break;
    pause;
end


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

%%
figure(1); set(gcf, 'Color', 'white', 'Position', [0 0 1000 300]);
h = boxplot(scores(:,:,3), 'orientation', 'horizontal');
set(h,{'linew'},{2});
xlabel('SD')
set(gca,'ytick',1:3);
set(gca,'yticklabel',{'neuTube', 'FullFlow', 'Ours'}, 'FontSize', 14);
xlim([0 30])
figure(2); set(gcf, 'Color', 'white', 'Position', [0 0 1000 300]);
h = boxplot(scores(:,:,5), 'orientation', 'horizontal');
set(h,{'linew'},{2});
xlabel('%SSD')
set(gca,'ytick',1:3);
set(gca,'yticklabel',{'neuTube', 'FullFlow', 'Ours'}, 'FontSize', 14);
