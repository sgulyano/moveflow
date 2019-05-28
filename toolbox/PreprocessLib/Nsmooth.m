function smoothed = Nsmooth(im,sigGauss)

%% smooth and find gradient
% im = medfilt2(im,[5 5]);
% hyfilt = fspecial('sobel');
% hxfilt = hyfilt';
% Iy1 = imfilter(double(im),hyfilt,'replicate');
% Ix1 = imfilter(double(im),hxfilt,'replicate');


% gauss_sigma = 1.0;
maskSize = ceil(sqrt(sigGauss*6));

    gmap = double(fspecial('gaussian', maskSize, sigGauss));
   
    gLeft  = [gmap, zeros(maskSize,1)];
    gRight = [zeros(maskSize,1), gmap];
    gDiff = gLeft-gRight;
    DoG.Gx = gDiff(:,1:maskSize);

    gTop    = [gmap; zeros(1,maskSize)];
    gBottom = [zeros(1,maskSize); gmap];
    gDiff = gTop-gBottom;
    DoG.Gy = gDiff(1:maskSize,:);
    
 Ix2 = conv2fft(double(im),DoG.Gx,'same');
 Iy2 = conv2fft(double(im),DoG.Gy,'same');

smoothed = sqrt(Ix2.^2+Iy2.^2);

%imshow(smoothed,[]);
