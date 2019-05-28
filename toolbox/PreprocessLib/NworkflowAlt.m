%function [curvelets,bkgrd,voting,sFlux]= Nworkflow(depthImg)
tStart=tic;
tic
[m,n,o]=size(depthImg);
voxels=m*n*o;
depthImg=single(depthImg);
thresh=depthImg<25;
depthImg(thresh)=0; %change all values below 25 to 0
clear thresh;
%% Curvelet Preprocessing
curvelets=zeros(size(depthImg));
sigCurv=20;
for i=1:o   
        disp(['Curvelets #',num2str(i),' started']);
        [curvelets(:,:,i)]=NfdctDemo(depthImg(:,:,i),sigCurv);
end
toc

%% Quick Salt-and-Pepper Cleanup with a Median Filter
tic
filtered = zeros(size(curvelets));

for i=1:o
    disp(['Median Filtering #',num2str(i),' started']);
    [filtered(:,:,i)]=medfilt2(curvelets(:,:,i),[3,3],'symmetric');
end
toc
clear curvelets;

%% Remove connected sets of voxels smaller than 200 elements
tic;
disp('Small Region Removal started');
filtered=filtered/max(filtered(:));
filtered=filtered.*255;
CC = bwconncomp(filtered);
numPixels = cellfun(@numel,CC.PixelIdxList);

for i=1:sum(numPixels<200)
    if numel(CC.PixelIdxList{i})<200
        filtered(CC.PixelIdxList{i})=0;
    end
end
clear CC;
toc

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

%% Block Frangi Vesselness Calc because the alternative is slow as %&^$.

blockTime=tic;
disp('Block Vesselness Calculation started');
options = struct('FrangiScaleRange',[1 3],'FrangiScaleRatio',2,'FrangiAlpha',.1,'FrangiBeta',10000,'FrangiC',100,'BlackWhite',false,'verbose',false);
Iout=zeros(size(filtered));
clear filtered;
for i=1:bRows
    for j=1:bCols %TODO Fix if input not evenly divisible by 32
        
        x=1+((j-1)*32);
        y=1+((i-1)*32);  
        [sIout,~,sL1,sL2,sL3,sVx,sVy,sVz]=NFrangiFilter3D(filteredPad(x:(x+31+2*pad),y:(y+31+2*pad),:),options,'n');
%     ,sVx2,sVy2,sVz2,sVx3,sVy3,sVz3
    Iout(x:(x+31),y:(y+31),:)=sIout((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    L1(x:(x+31),y:(y+31),:)=sL1((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    L2(x:(x+31),y:(y+31),:)=sL2((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    L3(x:(x+31),y:(y+31),:)=sL3((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    Vx(x:(x+31),y:(y+31),:)=sVx((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    Vy(x:(x+31),y:(y+31),:)=sVy((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    Vz(x:(x+31),y:(y+31),:)=sVz((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
%     Vx2(x:(x+31),y:(y+31),:)=sVx2((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
%     Vy2(x:(x+31),y:(y+31),:)=sVy2((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
%     Vz2(x:(x+31),y:(y+31),:)=sVz2((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
%     Vx3(x:(x+31),y:(y+31),:)=sVx3((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
%     Vy3(x:(x+31),y:(y+31),:)=sVy3((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
%     Vz3(x:(x+31),y:(y+31),:)=sVz3((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
    clear sIout sVx sVy sVz sVx2 sVy2 sVz2 sVx3 sVy3 sVz3
    
    end
    
end
clear filteredPad;
toc(blockTime);

%% Perform Block-Wise Graph Cuts
% Pare Down the Optimally Oriented Flux before Graph Cuts

Iout=Iout./max(Iout(:));
Iout=Iout.*255;
clear L1 L2 L3;
minCut=zeros(size(Iout));
blockTime=tic;
disp('Block Graph Cut started');
tau=20;
pad=1;
OOFOutPad=padarray(Iout,[pad pad pad],'symmetric');


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
Iout(~minCut)=0;
% clear minCut;

%% Find Seed Points Blockwise 
Iout=Iout./max(Iout(:));
Iout=Iout.*255;
tic; disp('Blockwise Seed Collection started');
fun=@(block_struct) sum(block_struct.data(:));
seedBlock=blockproc(Iout,[32,32],fun);
seeds=zeros(size(Iout));

for j=1:bCols
    for i=1:bRows %TODO Fix if input not evenly divisible by 32
        if seedBlock(i,j)> 0
            r=1+((i-1)*32);
            c=1+((j-1)*32); 
            block=Iout(r:(r+31),c:(c+31),:);
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
%%
clear m n o i j bCols bRows blockTime tStart


