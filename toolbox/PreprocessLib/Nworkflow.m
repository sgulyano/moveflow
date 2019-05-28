%function [curvelets,bkgrd,voting,sFlux]= Nworkflow(depthImg)

tStart=tic;

m=size(depthImg,1);
n=size(depthImg,2);
o=size(depthImg,3);
depthImg=double(depthImg);
curvelets=zeros(size(depthImg));
bkgrd=curvelets;
voting=bkgrd;

for i=1:o
    disp(['Curvelets #',num2str(i),' started']);
    tic, [curvelets(:,:,i)]=NfdctDemo(depthImg(:,:,i)); toc
    disp(['Background Subtraction #',num2str(i),' started']);
    tic, bkgrd(:,:,i)=Nsmooth(curvelets(:,:,i)); toc
    
end
bkgrdPad=cat(3,zeros(m,n),bkgrd,zeros(m,n));
disp('3D Spherical Flux Filtering');
tic, sFluxPad=Ngauss3filter(bkgrdPad); toc
sFlux=sFluxPad(:,:,2:24); 
disp('3D GVF Calculation');
tic, [Fext]=AM_GVF(sFlux,.1,20); toc
U=Fext(:,:,:,1);
V=Fext(:,:,:,2);
W=Fext(:,:,:,3);
disp('GVF Derivative Calculation');
tic
%%
tic
Dxx = Ngradient3(U,'x');
Dxy = Ngradient3(U,'y');
Dxz = Ngradient3(U,'z');
Dyy = Ngradient3(V,'y');
Dyz = Ngradient3(V,'z');
Dzz = Ngradient3(W,'z');
toc
%%

gradmag=sqrt(U.^2+V.^2+W.^2); %Include 3D data in what will end up being a 2D scalar voting
%%
cell=tic;
orientation = mod((V./U),pi);
for i=1:23
   disp(['Scalar Voting #',num2str(i),' started']);
   tic, voting(:,:,i)=NscalarVoting(2,1,Dxx(:,:,i),Dxy(:,:,i),Dyy(:,:,i),orientation(:,:,i),U(:,:,i),gradmag(:,:,i)); toc%,C,E);
end    
figure;
imshow3D(voting);
toc(cell)
toc(tStart)

