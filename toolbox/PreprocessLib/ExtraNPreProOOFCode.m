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

% clear bkgrd % cPad ~ curvelets with a 1 frame pad to avoid Fourier transform wrap-around artifacts
% tic; disp('3D Gaussian Smoothing started');
% smoothPad=Ngauss3filter(cPad); toc;
% clear bPad

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
% %% Optimally Oriented Flux from Preprocessed Image, 174s total
% 
% opts.useabsolute=1;
% opts.sigma=1;
% opts.responsetype=1; % 6 for vesselness eigenvalue calc, one for sum of two largest eigenvalues, i.e. where circular cross-sections exist
% opts.normalizationtype=1;
% disp(['Optimally Oriented Flux from preprocessed image started']);
% tic;   OOFOut=Noof3response(smoothPad,1:5,opts); toc
% OOFOut=OOFOut(:,:,3:o+2);

% disp(['Vesselness from GVF calculation started']);
% [Iout,Vx,Vy,Vz]=NFrangiGVF3D(ux,uy,uz,vy,vz,wz,options,m,n,o);
toc

% %% Perform 2D Graph Cuts
% disp('2D Graph Cut Calculation started');
% blockTime=tic;
% tau=20;
% %normalize the vesselness image before running graph cut
% OOFOut=OOFOut./max(OOFOut(:));
% OOFOut=OOFOut.*255;
% [minCut]=runGCO(OOFOut,tau,'2D');
% 
% toc(blockTime);
%% Perform a full 3D Graph Cut, a single block
% 
% tic
% disp('Full 3D Graph Cut Calculation started');
% tau=20;
% %normalize the vesselness image before running graph cut
% Iout=Iout./max(Iout(:));
% Iout=Iout.*255;
% minCut=runGCO(Iout,tau,'3D');
% imshow3D(minCut);
% toc
% runGCO(frangi,tau)
%% Centerline Calc

tic
disp(['Centerline calculation started']);
vec2dotU = Vx2.*fU;
vec2dotV = Vy2.*fV;
vec2dotW = Vz2.*fW;
vec3dotU = Vx3.*fU;
vec3dotV = Vy3.*fV;
vec3dotW = Vz3.*fW;
clear fU fV fW Vx2 Vy2 Vz2 Vx3 Vy3 Vz3
vec2dotGVF = vec2dotU+vec2dotV+vec2dotW;
vec3dotGVF = vec3dotU+vec3dotV+vec3dotW;
clear vec2dotU vec2dotV vec2dotW vec3dotU vec3dotV vec3dotW;
initialSeeds = ((vec2dotGVF+vec3dotGVF).^2)>.3;
clear vec2dotGVF vec3dotGVF;