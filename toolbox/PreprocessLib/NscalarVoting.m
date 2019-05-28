function voting=NscalarVoting(sigtensor,order,orientation,gradmag)%C,E)

%[Ix,Iy,orientation] = NgetCurvelet(im,sigcurvelets);%,C,E);
% % % gradmag = double(gradmag);
% % % Ix = 5*gradmag.*dirsin;
% % % Iy = -5*gradmag.*dircos;
% % 
% % %% compute the initial tensor
% % 
% Axx = Ix.*Ix ;%+ 0.001.*Iy.*Iy; % Ixp*Ixp
% Axy = Ix.*Iy ;%- 0.001.*Ix.*Iy;
% Ayy = Iy.*Iy ;%+ 0.001.*Ix.*Ix;
% 
% % orientation = mod(0.5*angle(Axx - Ayy + 2*1i*Axy)+pi,pi);

%Stickness should be the magnitude of the gradient, first-order style, so
%try 
stickness= gradmag;%sqrt(Dx.^2+Dy.^2); %gradmag=sqrt(Dx.^2+Dy.^2);

%stickness = sqrt((Axx+Ayy).^2 - 4*((Axx.*Ayy)-Axy.^2));
% % ballness = 0.5*((Axx+Ayy).^2 - stickness);
% 
%orientation = mod(atan2(dirsin,dircos),pi);%put in NgetCurvelet
% 
% % [N,X] = hist(stickness(:),100);
% % N = N./sum(N);
% % N = cumsum(N);
% %thresh = X(find(N>0.85,1));
% 
% % imagesc(orient*180/pi.*(stickness>0.5*10^4));
% % hold on;
% % [X,Y] = meshgrid(1:size(dircos,2),1:size(dircos,1));
% % quiver(X,Y,dircos.*stickness/max(stickness(:)),-dirsin.*stickness/max(stickness(:)));
% %  figure(1), imagesc(stickness);
% %  figure, imshow(uint8(128*ballness));
% %  figure(2), imagesc(orientation*180/pi.*(stickness>thresh));
% 
% % A0 = Axx + Ayy ;
% % A2 = Axx - 2*1i*Axy - Ayy;
% % A2n = Axx + 2*1i*Axy - Ayy;

sigma = sigtensor;
K = double(sigma*4);
[x,y] = meshgrid(-K:1:K,-K:1:K);
norm = sqrt(x.^2+y.^2);
norm(K+1,K+1) = 1;
mat = (x+1i*y)./norm;
mat = angle(mat);
mat = mod(mat,pi);
mat(K+1,K+1)=0;
% figure(4);
% order = 8;
a = [1 2 1];
b = [1 2 1];
for co = 1: order -1
    a = conv(a,b);
end
a = 2*a(1:(end+1)/2);
a(end) = a(end)/2;
a = a(end:-1:1);
for co = 0:2:4*order
    w{co/2+1} = exp(-(x.^2+y.^2)/2/sigma/sigma).*exp(-1i*mat*co);
%     w{co/2+1} = (x.^2+y.^2).*exp(1-(x.^2+y.^2)/2/sigma/sigma).*exp(-1i*mat*co);
%     subplot(5,2,(co/2)*2+1); imagesc(real(w{co/2+1}));
%     subplot(5,2,(co/2)*2+2); imagesc(imag(w{co/2+1}));
end
% figure(6);
for co = 0:2:4*order -2
    c{co/2+1} = stickness.*exp(-1i*co*orientation);
%     subplot(4,2,(co/2)*2+1); imagesc(real(c{co/2+1}));
%     subplot(4,2,(co/2)*2+2); imagesc(imag(c{co/2+1}));
end


if 1
% U2n = conv2fft(conj(c{2}),w{1},'same') + 4.0*conv2fft(c{1},w{2},'same')+6.0*conv2fft(c{2},w{3},'same') + 4*conv2fft(c{3},w{4},'same') + conv2fft(c{4},w{5},'same');
% U2 = conj(U2n);
U0 = real(a(1)*conv2fft(c{1},w{1},'same')); 
for co = 2:numel(a)
    U0 = U0 + real(a(co)*conv2fft(c{co},w{co},'same')); 
end
% U0 = real(6.0*conv2fft(c{1},w{1},'same') + 8.0*conv2fft(c{2},w{2},'same') + 2*conv2fft(c{3},w{3},'same'));
else
%     U2n = conv2(conj(c{2}),w{1},'same') + 4.0*conv2(c{1},w{2},'same')+6.0*conv2(c{2},w{3},'same') + 4*conv2(c{3},w{4},'same') + conv2(c{4},w{5},'same');
% U2 = conj(U2n);
U0 = real(6.0*conv2(c{1},w{1},'same') + 8.0*conv2(c{2},w{2},'same') + 2*conv2(c{3},w{3},'same'));
end
% sticknew = abs(U2n);
% orientnew = mod(0.5*angle(U2n)+pi,pi);
% ballnew = 0.5*(U0 - abs(U2));
% figure(3)
% subplot(2,2,1); imagesc(orientation*180/pi.*(stickness>thresh));
% subplot(2,2,2); imagesc(gradmag);
% subplot(2,2,3); imagesc(U0)
% subplot(2,2,4); imagesc(orientnew*180/pi.*(stickness>thresh));
% out = uint8(255*sticknew/max(sticknew(:)));
voting = U0;
% writeim(uint8(255*U0./max(U0(:))),'scalar_voting_output_3.tif');
% writeim(uint8(255*gradmag./max(gradmag(:))),'curvelet_output_3.tif');
end
% subplot(2,2,4); imagesc(ballnew);
