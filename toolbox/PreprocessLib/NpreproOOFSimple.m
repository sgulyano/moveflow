tStart=tic;
tic
[m,n,o]=size(depthImg);
voxels=m*n*o;
depthImg=double(depthImg);
curvelets=zeros(size(depthImg));
sigCurv=20;
for i=1:o
    
        disp(['Curvelets #',num2str(i),' started']);
        [curvelets(:,:,i)]=NfdctDemo(depthImg(:,:,i),sigCurv);
end
toc


% %% Adjusted Farsight Toolkit 2D Gaussian Smoothing
% sigGauss=80;
% bkgrd=zeros(size(depthImg));
% voting=bkgrd;
% 
% tic
% for i=1:o
%     disp(['Background Subtraction #',num2str(i),' started']);
%      bkgrd(:,:,i)=Nsmooth(curvelets(:,:,i),sigGauss);
% end
% toc
% % one=bkgrd(:,:,1,:);
% % one=cat(3,depthImg(:,:,1),one(:,:,1,1),one(:,:,1,2),one(:,:,1,3),one(:,:,1,4),one(:,:,1,5),one(:,:,1,6));
% % imshow3D(one);  %take a look at the first page of each of the different
% % calculated filtered images
% toc

%% 3D Gaussian smoothing using Fourier coefficients
% cPad=cat(3,zeros(m,n,2),curvelets,zeros(m,n,2));
% clear bkgrd % cPad ~ curvelets with a 1 frame pad to avoid Fourier transform wrap-around artifacts
% tic; disp('3D Gaussian Smoothing started');
% smoothPad=Ngauss3filter(cPad); toc;
% clear bPad
% smooth=cPad(:,:,3:o+2);

% %% Scalar Voting
% 
% 
% % blockRows=floor(m./32);
% % blockCols=floor(n./32);
% % blockDepth=o;
% 
% 
% % bkgrdCell=mat2cell(bkgrd,[blockRows],[blockCols],[blockDepth])
% % fun=@(bkgrd) AM_GVF(bkgrd,.1,20);
% %  [Fext]=blockproc3(depthImg,[blockRows blockCols blockDepth],fun);
%  %(Image,mu,iterations)
%  
% tic
% sigtensor=[1,2,4,8];
% order=1;
% for i=1:23
%     for j=1:4
% Dx = Ngradient3(bkgrd(:,:,i),'x');
% Dy = Ngradient3(bkgrd(:,:,i),'y');
% % Dxx = Ngradient3(Dx,'x');
% % Dxy = Ngradient3(Dx,'y');
% % Dyy = Ngradient3(Dy,'y');
% orientation= mod((Dx./Dy),pi);
% gradmag=sqrt(Dx.^2+Dy.^2);
% clear Dx Dy;
% disp(['Scalar Voting #',num2str(i),' started']);
% voting(:,:,i,j)=NscalarVoting(sigtensor(j),order,orientation,gradmag);%,C,E);
% clear orientation gradmag; %free memory for next round
%     end
% end
% toc

% %% Noise estimate from Frangi's S parameter
% tic
% disp('Noise from Frangi''s second-order structureness started');
% options = struct('FrangiScaleRange',[1,2],'FrangiScaleRatio',2,'BlackWhite',false);
% invNoise=NFrangiFilter3D(smooth,options,'s'); %NFrangiFilter3D(Image,a struct of options,type: any of 's','a','b','n', and n is the standard one)
% % Determine mu for each block based on amount of noise present in each
% % block.  invNoise is a measure of non-noise used in calculation of
% % vesselness
% % Block Processing of Gradient Vector Flow (GVF)
% 
% noise=(max(invNoise(:))+min(invNoise(:)))-invNoise;
% fun=@(block_struct) sum(block_struct.data(:));
% noiseBlock=blockproc(noise,[32,32],fun); %find level of noise in each 3D block
% numBlocks=numel(noiseBlock);
% [bRows,bCols]=size(noiseBlock);
% normalized=noiseBlock/max(noiseBlock(:));
% %Convert the normalized noise to a suitable range of mu: 0.0-0.2.  Mu is
% %the smoothness parameter for the GVF, and should be higher in noisier
% %areas.
% mus=normalized.*.2;
% toc
%%
% fU=zeros(size(smooth));
% fV=fU;
% fW=fV;
% tic
% disp('Block GVF started');
% for i=1:bRows
%     for j=1:bCols%TODO Fix if input not evenly divisible by 32
%         
%         x=1+((j-1)*32);
%         y=1+((i-1)*32);  
%     GVF=AM_GVF(smooth(x:(x+31),y:(y+31),:),mus(i,j),10);
%     fU(x:(x+31),y:(y+31),:)=GVF(:,:,:,1);
%     fV(x:(x+31),y:(y+31),:)=GVF(:,:,:,2);
%     fW(x:(x+31),y:(y+31),:)=GVF(:,:,:,3);
%     clear GVF U V W
%     end
% end
% toc
% %%
% 
% % ux=Ngradient3(fU,'x');
% % uy=Ngradient3(fU,'y');
% % uz=Ngradient3(fU,'z');
% % clear fU;
% % 
% % %vx=Ngradient3(V,'x');
% % vy=Ngradient3(fV,'y');
% % vz=Ngradient3(fV,'z');
% % clear fV;
% % %wx=Ngradient3(W,'x');
% % %wy=Ngradient3(W,'y');
% % wz=Ngradient3(fW,'z');
% % clear fW;
% % %%
% % % tic
% % % ux=reshape(ux,1,1,voxels);
% % % uy=reshape(uy,1,1,voxels);
% % % uz=reshape(uz,1,1,voxels);
% % % %vx=reshape(vx,1,1,voxels);
% % % vy=reshape(vy,1,1,voxels);
% % % vz=reshape(vz,1,1,voxels);
% % % %wx=reshape(wx,1,1,voxels);
% % % %wy=reshape(wy,1,1,voxels);
% % % wz=reshape(wz,1,1,voxels);
% % % jacobian=[ux,uy,uz;vx,vy,vz;wx,wy,wz];
% % % toc
% % 
% % % %% Eigenvalue calculation
% % % % tic; eigs=eig3(jacobian); toc
% % % %%
% % % tic;
% % % for i=1:length(jacobian)
% % %     eigsM(:,:,i)=eig(jacobian(:,:,i));
% % % end
% % % toc
% %% Optimally Oriented Flux from Preprocessed Image, 174s total
% 
% opts.useabsolute=1;
% opts.sigma=1;
% opts.responsetype=1; % 6 for vesselness eigenvalue calc, one for sum of two largest eigenvalues, i.e. where circular cross-sections exist
% opts.normalizationtype=1;
% disp(['Optimally Oriented Flux from preprocessed image started']);
% tic;   OOFOut=Noof3response(smoothPad,1:5,opts); toc
% OOFOut=OOFOut(:,:,3:o+2);

%% Alternatively Block-wise OOF (compared to full, only .00048% of the pixels are different (115 of 1024*1024*23) 77s total

blockTime=tic;
OOFOut=zeros(size(curvelets));
disp('Block Optimally Oriented Flux started');
for i=1:bRows
    for j=1:bCols%TODO Fix if input not evenly divisible by 32
        x=1+((j-1)*32);
        y=1+((i-1)*32);  
    [output,sL1,sL2,sL3,sVx,sVy,sVz,sVx2,sVy2,sVz2,sVx3,sVy3,sVz3]=Noof3response(curvelets(x:(x+31),y:(y+31),:),1:2,opts);
    OOFOut(x:(x+31),y:(y+31),:)=output;
    L1(x:(x+31),y:(y+31),:)=sL3;
    L2(x:(x+31),y:(y+31),:)=sL2;
    L3(x:(x+31),y:(y+31),:)=sL3;
%     Vx(x:(x+31),y:(y+31),:)=sVx;
%     Vy(x:(x+31),y:(y+31),:)=sVy;
%     Vz(x:(x+31),y:(y+31),:)=sVz;
%     Vx2(x:(x+31),y:(y+31),:)=sVx2;
%     Vy2(x:(x+31),y:(y+31),:)=sVy2;
%     Vz2(x:(x+31),y:(y+31),:)=sVz2;
%     Vx3(x:(x+31),y:(y+31),:)=sVx3;
%     Vy3(x:(x+31),y:(y+31),:)=sVy3;
%     Vz3(x:(x+31),y:(y+31),:)=sVz3;
    clear output sL1 sL2 sL3 sVx sVy sVz sVx2 sVy2 sVz2 sVx3 sVy3 sVz3
   end
end
% OOFOut=OOFOut(:,:,3:o+2);
% Vx=Vx(:,:,3:o+2);
% Vy=Vy(:,:,3:o+2);
% Vz=Vz(:,:,3:o+2);
% Vx2=Vx2(:,:,3:o+2);
% Vy2=Vy2(:,:,3:o+2);
% Vz2=Vz2(:,:,3:o+2);
% Vx3=Vx3(:,:,3:o+2);
% Vy3=Vy3(:,:,3:o+2);
% Vz3=Vz3(:,:,3:o+2);
toc(blockTime);

%% Centerline Calc
% clear Vx Vy Vz
% tic
% disp(['Centerline calculation started']);
% vec2dotU = Vx2.*fU;
% vec2dotV = Vy2.*fV;
% vec2dotW = Vz2.*fW;
% vec3dotU = Vx3.*fU;
% vec3dotV = Vy3.*fV;
% vec3dotW = Vz3.*fW;
% clear fU fV fW Vx2 Vy2 Vz2 Vx3 Vy3 Vz3
% vec2dotGVF = vec2dotU+vec2dotV+vec2dotW;
% vec3dotGVF = vec3dotU+vec3dotV+vec3dotW;
% clear vec2dotU vec2dotV vec2dotW vec3dotU vec3dotV vec3dotW;
% centerlines = ((vec2dotGVF+vec3dotGVF).^2)<.00001;
% clear vec2dotGVF vec3dotGVF;
% % disp(['Vesselness from GVF calculation started']);
% % [Iout,Vx,Vy,Vz]=NFrangiGVF3D(ux,uy,uz,vy,vz,wz,options,m,n,o);
% toc

%% Perform 2D Graph Cuts
disp('2D Graph Cut Calculation started');
blockTime=tic;
tau=20;
%normalize the vesselness image before running graph cut
OOFOut=OOFOut./max(OOFOut(:));
OOFOut=OOFOut.*255;
[minCut]=runGCO(OOFOut,tau,'2D');

toc(blockTime);

% %% Perform Block-Wise Graph Cuts
% minCut=zeros(size(depthImg));
% blockTime=tic;
% disp('Block Graph Cut Calculation started');
% tau=20;
% %normalize the vesselness image before running graph cut
% 
% OOFOut=OOFOut./max(OOFOut(:));
% OOFOut=OOFOut.*255;
% for i=1:bRows
%     for j=1:bCols %TODO Fix if input not evenly divisible by 32
%         x=1+((j-1)*32);
%         y=1+((i-1)*32);  
%         miniCut=runGCO(OOFOut(x:(x+31),y:(y+31),:),tau,'3D');
%         minCut(x:(x+31),y:(y+31),:)=miniCut;
%         clear miniCut;
%         clear gco_matlab;
%     end
% end
% % clear OOFOut x y
% toc(blockTime);


%% Find Seed Points
seeds=minCut&centerlines;
% clear minCut centerlines;
% imshow3D(seeds);
 
toc(tStart)
clear m n o i j bCols bRows blockTime tStart