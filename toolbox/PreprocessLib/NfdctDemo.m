% %disp(' ');
% %disp('fdct_wrapping_demo_denoise.m -- Image denoising using Curvelets');
% %disp(' ');
% %disp('Denoising is achieved by hard-thresholding of the curvelet coefficients.');
% %disp('We select the thresholding at 3*sigma_jl for all but the finest scale');
% %disp('where it is set at 4*sigmajl; here sigma_jl is the noise level of a');
% %disp('coefficient at scale j and angle l (equal to the noise level times');
% %disp('the l2 norm of the corresponding curvelet). There are many ways to compute');
% %disp('the sigma_jl''s, e.g. by computing the norm of each individual curvelet,');
% %disp('and in this demo, we do an exact computation by applying a forward curvelet');
% %disp('transform on an image containing a delta function at its center.');
% %disp(' ');

% fdct_wrapping_demo_denoise.m -- Image denoising using Curvelets
function [fdctImg]=NfdctDemo(img,sigma)
img=double(img);
m=size(img,1);
n=size(img,2);
% img = double(imread('Lena.jpg'));
%scales=ceil(log2(min(m,n)) - 4);

% sigma = 20;        
is_real = 0;

noisy_img = img; %+ sigma*randn(n);

%disp('Compute all thresholds');
F = ones(m,n);
X = fftshift(ifft2(F)) * sqrt(numel(F));
 C = NfdctWrapping(X,is_real,2); 
%curveletC=C;
% Compute norm of curvelets (exact)
E = cell(size(C));
for s=1:length(C)
  E{s} = cell(size(C{s}));
  for w=1:length(C{s})
    A = C{s}{w};
    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
  end
end
%curveletE=E;
% Take curvelet transform
%disp(' ');
%disp('Take curvelet transform: fdct_wrapping');
 C = NfdctWrapping(noisy_img,is_real,2); 

% Apply thresholding
Ct = C;
for s = 2:length(C)
  thresh = 3*sigma + sigma*(s == length(C));
  for w = 1:length(C{s})
    Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
  end
end

% Take inverse curvelet transform 
%disp(' ');
%disp('Take inverse transform of thresholded data: ifdct_wrapping');
 fdctImg = real(NifdctWrapping(Ct,1,m,n)); 
%figure;
%imshow(fdctImg,[]);

% subplot(1,3,1); imagesc(img); colormap gray; axis('image');
% subplot(1,3,2); imagesc(noisy_img); colormap gray; axis('image');
% subplot(1,3,3); imagesc(restored_img); colormap gray; axis('image');
% %%
% %disp(' ');
% %disp('This method is of course a little naive.');
% R=input('Do you want to see the outcome of a more sophisticated experiment? [Y/N] ','s');
% %disp(' ');
% 
% if strcmp(R,'') + strcmp(R,'y') + strcmp(R,'Y'), 
%   combined = double(imread('LenaCombined.jpg'));
%   figure; imagesc(combined); colormap gray; axis('image');
%   
% %disp('This image (courtesy of Jean-Luc Starck) is the result of a more sophisticated');
% %disp('strategy which involves both curvelets and wavelets. For a reference, please see');
% %disp('J.L. Starck, D.L. Donoho and E. Candes, Very High Quality Image Restoration,')
% %disp('in SPIE conference on Signal and Image Processing: Wavelet Applications');
% %disp('in Signal and Image Processing IX, A. Laine, M. A. Unser and A. Aldroubi Eds,'); 
% %disp('Vol 4478, 2001.');
% %disp(' ');
% %disp('For other experiments, please check fdct_wrapping_demo_denoise_enhanced.');
% %disp(' ')
% end
