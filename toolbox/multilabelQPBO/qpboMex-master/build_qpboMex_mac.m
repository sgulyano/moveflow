function build_qpboMex_mac
% build_qpboMex builds package qpboMex for Mac based on solution from
% http://stackoverflow.com/questions/15298603/linking-errors-from-matlab-mex-library
%
% Sarun Gulyanon,  07.04.2017

directory = mfilename('fullpath'); directory(end-16:end) = [];
codePath = [directory 'QPBO-v1.32.src'];
srcFiles = 'qpboMex_mac.cpp';

cmdLine = ['mex ' directory srcFiles ' -output ' directory 'qpboMex -largeArrayDims ', '-I', codePath];
eval(cmdLine);




