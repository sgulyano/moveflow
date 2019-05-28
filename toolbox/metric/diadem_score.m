function [score, miss, extra] = diadem_score( goldfile, testfile, neurontype )
%EV_DIADEM find diadem score of trace in testfile compared to gold standard
%gold file

score = 0;
miss = zeros(0,4);
extra = zeros(0,4);

functionname='diadem_score.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
currentdir = pwd;
cd(functiondir);

bashfile = 'diadem.sh';
resultfile = 'result.txt';
if nargin < 3
    neurontype = 5;
end

if ~ischar(goldfile)
    goldswc = goldfile;
    goldfile = 'gold.swc';
    fileID = fopen(goldfile, 'w');
    fprintf(fileID,'%d %d %.3f %.3f %.3f %.4f %d\n',goldswc');
    fclose(fileID);
end
if ~ischar(testfile)
    testswc = testfile;
    testfile = 'test.swc';
    fileID = fopen(testfile, 'w');
    fprintf(fileID,'%d %d %.3f %.3f %.3f %.4f %d\n',testswc');
    fclose(fileID);
end

%% run DiademMetric
fileID = fopen(bashfile, 'w');
% this line to fix java problem in Mac/MATLAB
fprintf(fileID,'unset DYLD_FRAMEWORK_PATH DYLD_LIBRARY_PATH\n');
out = sprintf('java -jar DiademMetric.jar -G %s -T %s -D %d -m true > %s', ...
        goldfile, testfile, neurontype, resultfile);
fprintf(fileID,'%s\n', out);
fclose(fileID);

system(['bash ' bashfile]);

%% read and parse result
fid = fopen(resultfile, 'r');
tline = fgets(fid);
while ischar(tline)
    % get diadem metric score
    tmp = regexp(tline,'Score:.(.*)','tokens');
    if ~isempty(tmp)
        score = str2double([tmp{1}{:}]);
    end
    
    % get missed nodes
    if strncmp(tline,'Nodes that were missed',22)
        cur = 1;
        tline = fgets(fid);
        if ischar(tline)
            tmp = regexp(tline,'(\d+\.\d+|\d+)','match');
            while ~isempty(tmp)
                miss(cur,:) = cellfun(@str2double, tmp);
                tline = fgets(fid);
                if ~ischar(tline), break, end
                tmp = regexp(tline,'(\d+\.\d+|\d+)','match');
                cur = cur + 1;
            end
        end
    end
    
    % get extra nodes
    if strncmp(tline,'Extra nodes in test reconstruction',34)
        cur = 1;
        tline = fgets(fid);
        if ischar(tline)
            tmp = regexp(tline,'(\d+\.\d+|\d+)','match');
            while ~isempty(tmp)
                extra(cur,:) = cellfun(@str2double, tmp);
                tline = fgets(fid);
                if ~ischar(tline), break, end
                tmp = regexp(tline,'(\d+\.\d+|\d+)','match');
                cur = cur + 1;
            end
        end
    end

    tline = fgets(fid);
end
fclose(fid);
delete(bashfile, resultfile);
if exist('goldswc', 'var'), delete(goldfile); end;
if exist('testswc', 'var'), delete(testfile); end;
cd(currentdir)
end

