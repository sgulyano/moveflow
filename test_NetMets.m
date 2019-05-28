close all; clc;

% filepath where vaa3d is installed
netmets = '~/Desktop/netmets.git/build/netmets';
directory = '~/Desktop/LipingsData/';
switch(1)
    case 1
        tracefile = 'swc/%s/larva3_%s_n%d_t%03d.swc';
        golddir = 'GoldStandard/ZstackL1_3_2grayscale/';  % directory containing the gold standard
        initswctime = [198, 201, 229, 252, 340, 342];
        outputfile = 'larva3';
        xlsname = 'larva3_NetMetsdata.xlsx';
    case 2
        tracefile = 'swc/%s/larva4S14Z_%s_n%d_t%03d.swc';
        golddir = 'GoldStandard/larva4S14Z/';
        initswctime = 1;
        outputfile = 'larva4S14Z';
        xlsname = 'larva4S14Z_NetMetsdata.xlsx';
end
methods = {'neutubeone', 'fullflow', 'our'};

if ~exist([outputfile '.txt'], 'file')
    %% read files
    swcname = dir([directory golddir '*.swc']);
    swcname(arrayfun(@(x)(x.name(1)=='.'), swcname)) = [];

    neuron_num = arrayfun(@(x)regexp(x.name,'.*_n(\d+)_t\d+.swc','tokens'), swcname);
    neuron_num = cellfun(@(x)str2double(x{1}), neuron_num);

    swc_time = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swcname);
    swc_time = cellfun(@(x)str2double(x{1}), swc_time);

    %% compute netmets
    bashfile = [outputfile '.sh'];
    for i = 1:length(neuron_num)
        num = neuron_num(i);
        t = swc_time(i);
        if t == initswctime(num), continue, end;
        goldfile = [golddir swcname(i).name];
        goldswc = readswc([directory golddir swcname(i).name]);
        score = zeros(length(methods), 2);
        for j = 1:length(methods)
            testfile = sprintf(tracefile, methods{j}, methods{j}, num, t);
            testswc = readswc(testfile);
            if isempty(testswc)
                testfile = 'swc/swcempty.swc'; 
            end

            if i == 1 && j == 1
                fileID = fopen(bashfile, 'w');
                fprintf(fileID,'echo %s_%d_%d > %s\n', methods{j}, num, t, [outputfile '.txt']);
            else
                fileID = fopen(bashfile, 'a');
                fprintf(fileID,'echo %s_%d_%d >> %s\n', methods{j}, num, t, [outputfile '.txt']);
            end
            fprintf(fileID,'%s %s %s --sigma 4 >> %s\n', netmets, goldfile, testfile, [outputfile '.txt']);
            fclose(fileID);
        end
    end
    disp('Done');
else
    fn = [];
    fp = [];
    nums = [];
    ts = [];

    fileID = fopen([outputfile '.txt'], 'r');
    tline = fgetl(fileID);
    while ischar(tline)
        tokens=strsplit(tline,'_');
        method = tokens{1};
        num = str2double(tokens{2});
        t = str2double(tokens{3});
        
        tline = fgetl(fileID);
        tokens = strsplit(tline, ' ');
        FNR = str2double(tokens{2});
        
        tline = fgetl(fileID);
        tokens = strsplit(tline, ' ');
        FPR = str2double(tokens{2});
        
        nums = [nums, num];
        ts = [ts, t];
        fn = [fn, FNR];
        fp = [fp, FPR];
        fprintf('%15s num=%d, t=%d => FNR: %f, FPR: %f\n', method, num, t, FNR, FPR);
        tline = fgetl(fileID);
    end

    nums = nums(1:3:end);
    ts = ts(1:3:end);
    fp = reshape(fp, 3, length(fp)/3);
    fn = reshape(fn, 3, length(fn)/3);

    T = array2table([nums; ts; fn]');
    T.Properties.VariableNames = [{'Frame'},{'NeuronNumber'},methods];
    writetable(T, xlsname,'Sheet','FNR');
    T = array2table([nums; ts; fp]');
    T.Properties.VariableNames = [{'Frame'},{'NeuronNumber'},methods];
    writetable(T, xlsname,'Sheet','FPR');
    %%
    figure(1); set(gcf, 'Color', 'w', 'Position', [0 0 1000 300]);
    h = boxplot(fn', 'Orientation', 'horizontal');
    set(h,{'linew'},{2})
    xlabel('FNR');
    set(gca,'ytick',1:3);
    set(gca,'yticklabel',{'neuTube', 'FullFlow', 'Ours'}, 'FontSize', 14);
    figure(2); set(gcf, 'Color', 'w', 'Position', [0 0 1000 300]);
    h = boxplot(fp', 'Orientation', 'horizontal');
    set(h,{'linew'},{2})
    xlabel('FPR');
    set(gca,'ytick',1:3);
    set(gca,'yticklabel',{'neuTube', 'FullFlow', 'Ours'}, 'FontSize', 14);
end
