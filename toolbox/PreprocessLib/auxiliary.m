% % function automatedLoadRunShow(depthImg)
% i=3;
% j=6;
% for a = [.05,.15,.25,.35,.45]
%     i=i+1
%     j=j-6
%     for c = [70 100 130 160 190]
%         j=j+1
%         options = struct('FrangiScaleRange',[1,10],'FrangiScaleRatio',1,'FrangiAlpha',a,'FrangiBeta',10000,'FrangiC',c,'BlackWhite',false,'verbose',true)
%         [Iout,whatScale] = FrangiFilter3D(depthImg,options);
%         fname = strcat('Iout',num2str(i),num2str(j),'.mat')
%         save(fname, 'Iout','whatScale');
%         clear Iout whatScale
%     end
% end
% %% Find Seed Points Blockwise
% OOFOut=OOFOut./max(OOFOut(:));
% OOFOut=OOFOut.*255;
% tic; disp('Blockwise Seed Collection started');
% fun=@(block_struct) sum(block_struct.data(:));
% seedBlock=blockproc(OOFOut,[32,32],fun);
% seeds=zeros(size(OOFOut));
% for j=1:bCols
%     for i=1:bRows %TODO Fix if input not evenly divisible by 32
%         
%         if seedBlock(i,j)> 0
%             x=1+((i-1)*32);
%             y=1+((j-1)*32); 
%             block=OOFOut(x:(x+31),y:(y+31),:);
%             miniSeeds=imregionalmax(block);
%             seeds(x:(x+31),y:(y+31),:)=miniSeeds;
%             clear miniSeeds;
% %         else
% %             seeds(x:(x+31),y:(y+31),:)=zeros(32,32,o);
% %             clear miniSeeds;
%         end
%     end
% end
% seedBlock=blockproc(seeds,[32,32],fun);
% toc
% response=L2+L3;
% response(L2<=0)=0;
% response(L3<=0)=0;
% response(response<0)=0;
% response=response/max(response(:));
% response=response.*255;
% 
% CC = bwconncomp(response);
% numPixels = cellfun(@numel,CC.PixelIdxList);
% 
% for i=1:sum(numPixels<200)
%     if numel(CC.PixelIdxList{i})<200
%         response(CC.PixelIdxList{i})=0;
%     end
% end
% response=response.*(1-abs(L1)/abs(L2));
% CC = bwlabeln(tFilt,6);
% % numPixels = cellfun(@numel,CC.PixelIdxList);
% % [biggest,idx] = max(numPixels);
% % largest=(CC.PixelIdxList{idx});
% % logicalL=zeros(size(depthImg));
% % logicalL(largest)=tFilt(largest);
% % imshow3D(logicalL);
% 
% response=response/max(response(:));
% response=response.*255;
% thresh=response<25;
% tFilt=response;
% tFilt(thresh)=0;


% for i=1:sum(numPixels<200)
%     if numel(CC.PixelIdxList{i})<200
%         tFilt(CC.PixelIdxList{i})=0;
%     end
% end
% imshow3D(tFilt);

%  minifU=fU(382:581,1:200,:);
%  minifV=fV(382:581,1:200,:);
%  minifW=fW(382:581,1:200,:);
%  miniVx=Vx(382:581,1:200,:);
%  miniVx2=Vx2(382:581,1:200,:);
%  miniVx3=Vx3(382:581,1:200,:);
%  miniVy=Vy(382:581,1:200,:);
%  miniVy2=Vy2(382:581,1:200,:);
%  miniVy3=Vy3(382:581,1:200,:);
%  miniVz=Vz(382:581,1:200,:);
%  miniVz2=Vz2(382:581,1:200,:);
%  miniVz2=Vz2(382:581,1:200,:);
%  miniL1=L1(382:581,1:200,:);
%  miniL2=L2(382:581,1:200,:);
%  miniL3=L3(382:581,1:200,:);
%  miniMinCut=minCut(382:581,1:200,:);  

% for i=1:23
%     i
%     OOFFiltered(:,:,i)=medfilt2(OOFOut(:,:,i),[4,4],'symmetric');
% end

% for i = 1:6
%     extra=i/6
%     output(:,:,:,i)=Noof3response(mini,1:.5:5,opts,extra);
% end


% [Fext]=AM_GVF(result2,.1,20);
% U=Fext(:,:,:,1);
% V=Fext(:,:,:,2);
% W=Fext(:,:,:,3);
% Dxx = Ngradient3(U,x);
% Dxy = Ngradient3(U,y);
% Dxz = Ngradient3(U,z);
% Dyy = Ngradient3(V,y);
% Dyz = Ngradient3(V,z);
% Dzz = Ngradient3(W,z);
% orientation = tan(V/U);
% invim=(max(Iouts(:))+min(Iouts(:)))-Iouts;
% I=Ioutn;
% fun=@(block_struct) sum(block_struct.data(:));
% I2D=sum(I,3);
% noise=blockproc(I2D,[32,32],fun);
% imshow(noise,[]);
% for i =1:23
%     bw(:,:,i)=im2bw(response(:,:,i),mean(response(:)));
% end
% imshow3D(bw);

% opts.useabsolute=1;
% opts.sigma=1;
% tStart=tic;
% for i =1:3
%     opts.normalizationtype=(i-1);
%    
%     for j=1:6
%          opts.responsetype=(j-1);
%       tic;   OOFOut=oof3response(curvelets,1:5,opts); toc
%         fname = strcat('OOFOut',num2str(i),num2str(j),'.mat')
%         save(fname, 'OOFOut');
%         clear OOFOut
%     end
% end
% toc(tStart);
% count=0;
% for i =1:3
%     for j=1:6
%         count=count+1
%         fname = strcat('OOFOut',num2str(i),num2str(j),'.mat');
%         load(fname,'OOFOut');
%         slices(:,:,count)=OOFOut(:,:,10);
%     end
% end
% 
% for i =1:23
% =quiver(vx3(:,:,i),vy3(:,:,i));
% end
% imshow3D(h);