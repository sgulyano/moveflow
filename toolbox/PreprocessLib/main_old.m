% function [seeds,minCut,OOFOut]=NpreproOOF(depthImg)
% author Nates

clear all
close all
clc

addpath(genpath('./CurveLab-2.1.3'))
addpath(genpath('./gco-v3.0'))

set(0,'DefaultFigureWindowStyle','docked')

%% Read Image
disp('Reading images ...')
currentpath = cd;
[filename, datapath] = uigetfiles;
cd(datapath)
im_names = dir('*.tif');
FNUM = length(im_names);
for fcount = 1:FNUM,
    X = imread(im_names(fcount).name,'tif');
    depthImg(:,:,fcount)=X;
end

% get name of save file
C = strsplit(deblank(datapath),'/');
C = C(~cellfun(@(x)isempty(x), C));
savefile = ['n', C{end}, '.mat'];

Vimg = depthImg;
cd(currentpath)
figure(1)
imshow3D(depthImg)
%%
tStart=tic;
tic
[m,n,o]=size(depthImg);
voxels=m*n*o;
depthImg=single(depthImg);
thresh=depthImg<25;
depthImg(thresh)=0; %change all values below 25 to 0
figure(2)
imshow3D(thresh)
clear thresh;
%% Curvelet Preprocessing
curvelets=zeros(size(depthImg));
sigCurv=20;
for i=1:o   
        disp(['Curvelets #',num2str(i),' started']);
        [curvelets(:,:,i)]=NfdctDemo(depthImg(:,:,i),sigCurv);
end
toc
figure(3)
imshow3D(curvelets)
%% Quick Salt-and-Pepper Cleanup with a Median Filter
tic
filtered = zeros(size(curvelets));

for i=1:o
    disp(['Median Filtering #',num2str(i),' started']);
    [filtered(:,:,i)]=medfilt2(curvelets(:,:,i),[3,3],'symmetric');
end
toc
clear curvelets;
figure(4)
imshow3D(filtered)
%% Remove connected sets of voxels smaller than 200 elements
tic;
disp('Small Region Removal started');
filtered=filtered/max(filtered(:));
filtered=filtered.*255;
CC = bwconncomp(filtered);
numPixels = cellfun(@numel,CC.PixelIdxList);

for i=1:length(numPixels)
    if numel(CC.PixelIdxList{i})<200
        filtered(CC.PixelIdxList{i})=0;
    end
end
clear CC;
toc
figure(5)
imshow3D(filtered)
%% Block GVF
fun=@(block_struct) sum(block_struct.data(:));
noiseBlock=blockproc(filtered,[32,32],fun);
[bRows,bCols]=size(noiseBlock);
fU=zeros(size(filtered));
fV=fU;
fW=fV;
pad=1;
filteredPad=padarray(filtered,[pad pad pad],'symmetric');

tic
disp('Block GVF started');

for j=1:bCols
    for i=1:bRows%TODO Fix if input not evenly divisible by 32  
        r=1+((i-1)*32);
        c=1+((j-1)*32);  
        block=filteredPad(r:(r+31+2*pad),c:(c+31+2*pad),:);
        [GVF,energy]=AM_GVF_en(block,.1,10);
        fU(r:(r+31),c:(c+31),:)=GVF((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad),1);
        fV(r:(r+31),c:(c+31),:)=GVF((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad),2);
        fW(r:(r+31),c:(c+31),:)=GVF((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad),3);
        fEnergy(r:(r+31),c:(c+31),:)=energy((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
        clear GVF U V W energy
    end
end
toc
clear filteredPad

%% Block-wise OOF (compared to full, only .00048% of the pixels are different (115 of 1024*1024*23) 77s total
blockTime=tic;
disp('Block Optimally Oriented Flux started');
OOFOut=zeros(size(filtered));
opts.useabsolute=1;
opts.sigma=1;
opts.responsetype=1; % 6 for vesselness eigenvalue calc, one for sum of two largest eigenvalues, i.e. where circular cross-sections exist
opts.normalizationtype=1;
pad=8;
opts.marginwidth=[pad pad pad];
filteredPad=padarray(filtered,opts.marginwidth,'symmetric');
% opts=rmfield(opts,'marginwidth');
% clear filtered;
for j=1:bCols
    for i=1:bRows%TODO Fix if input not evenly divisible by 32
        r=1+((i-1)*32);
        c=1+((j-1)*32);
        block=single(filteredPad(r:(r+31+2*opts.marginwidth(2)),c:(c+31+2*opts.marginwidth(1)),:));
        [output,sL1,sL2,sL3,sVx,sVy,sVz]=Noof3response(block,1:.5:7,opts);
        OOFOut(r:(r+31),c:(c+31),:)=output;
        L1(r:(r+31),c:(c+31),:)=sL1;
        L2(r:(r+31),c:(c+31),:)=sL2;
        L3(r:(r+31),c:(c+31),:)=sL3;
        Vx(r:(r+31),c:(c+31),:)=sVx;
        Vy(r:(r+31),c:(c+31),:)=sVy;
        Vz(r:(r+31),c:(c+31),:)=sVz;
        clear output sL1 sL2 sL3 sVx sVy sVz sVx2 sVy2 sVz2 sVx3 sVy3 sVz3
    end
end
toc(blockTime);
clear filteredPad


%% Perform Block-Wise Graph Cuts
% Pare Down the Optimally Oriented Flux before Graph Cuts
OOFOut(L2<=0)=0;
OOFOut(L3<=0)=0;
OOFOut(OOFOut<0)=0;
% OOFOut=OOFOut.*(1-abs(L1)./abs(L2));%scale according to ratio from Sato et al L1/L2 that should be closer to 0 for line-like structures
OOFOut=OOFOut./max(OOFOut(:));
OOFOut=OOFOut.*255;
clear L1 L2 L3;
minCut=zeros(size(OOFOut));
blockTime=tic;
disp('Block Graph Cut started');
tau=20;
pad=1;
OOFOutPad=padarray(OOFOut,[pad pad pad],'symmetric');
%normalize the oriented flux image before running graph cut
% OOFOutPad=OOFOutPad./max(OOFOutPad(:));
% OOFOutPad=OOFOutPad.*255;

for j=1:bCols
    for i=1:bRows %TODO Fix if input not evenly divisible by 32
        r=1+((i-1)*32);
        c=1+((j-1)*32);  
        block=OOFOutPad(r:(r+31+2*pad),c:(c+31+2*pad),:);
        miniCut=runGCO(block,tau,'3D');
        minCut(r:(r+31),c:(c+31),:)=miniCut((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
        clear miniCut;
        clear gco_matlab;
    end
end
clear OOFOutPad
toc(blockTime);
minCut=logical(minCut);
OOFOut(~minCut)=0;
% clear minCut;

%% Find Seed Points Blockwise 
OOFOut=OOFOut./max(OOFOut(:));
OOFOut=OOFOut.*255;
tic; disp('Blockwise Seed Collection started');
fun=@(block_struct) sum(block_struct.data(:));
seedBlock=blockproc(OOFOut,[32,32],fun);
seeds=zeros(size(OOFOut));

for j=1:bCols
    for i=1:bRows %TODO Fix if input not evenly divisible by 32
        if seedBlock(i,j)> 0
            r=1+((i-1)*32);
            c=1+((j-1)*32); 
            block=OOFOut(r:(r+31),c:(c+31),:);
            miniSeeds=imregionalmax(block);
            seeds(r:(r+31),c:(c+31),:)=miniSeeds;
            clear miniSeeds;
%       else
%           seeds(r:(r+31),c:(c+31),:)=zeros(32,32,o);
%           clear miniSeeds;
        end
    end
end
seedBlock=blockproc(seeds,[32,32],fun);
toc
% clear minCut centerlines;
% imshow3D(seeds);

toc(tStart)

%% Save filtered image and seed points
seeds = logical(seeds);
smooth = single(filtered);
V = Vimg;
blksize = [32 32 FNUM];
save (savefile, 'V', 'smooth', 'seeds', 'blksize', '-v7.3');
disp('Done');

%% Clear unused varaibles
clear m n o i j bCols bRows blockTime tStart

%% Show seed points
figure(6)
for fcount = 1:FNUM
    fcount
    imshow(filtered(:,:,fcount),[])
    [py, px] = find(seeds(:,:,fcount));
    hold on;
    plot(px,py,'b+')
    hold off;
    pause(1)
end
