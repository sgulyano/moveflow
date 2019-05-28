function build_qpboMex
% build_qpboMex builds package qpboMex
%
% Anton Osokin (firstname.lastname@gmail.com),  24.09.2014
% Update by Sarun Gulyanon, 07.04.2017

directory = mfilename('fullpath'); directory(end-12:end) = [];
codePath = [directory 'QPBO-v1.32.src'];

srcFiles = {fullfile(directory, 'qpboMex.cpp'), ...
            fullfile(codePath, 'QPBO.cpp'), ...
            fullfile(codePath, 'QPBO_maxflow.cpp'), ...
            fullfile(codePath, 'QPBO_postprocessing.cpp'), ...
            fullfile(codePath, 'QPBO_extra.cpp') };
allFiles = '';
for iFile = 1 : length(srcFiles)
    allFiles = [allFiles, ' ', srcFiles{iFile}];
end

cmdLine = ['mex ', allFiles, ' -output qpboMex -largeArrayDims ', '-I', codePath];
eval(cmdLine);




