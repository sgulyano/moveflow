function [T, scale] = IV_frangi( V )
%IV_FRANGI Summary of this function goes here
%   Detailed explanation goes here
%% Compute Frangi to find ev
r = 8;

% opt.FrangiScaleRange = [1,r];
% opt.FrangiScaleRatio = 1;
% opt.FrangiAlpha = .5;
% opt.FrangiBeta = .5;
% opt.FrangiC = 500;
% opt.BlackWhite = false;
% opt.verbose = true;
% 
% % pad image
% Vnew = min(255, double(V) * 8);
% [T,scale] = FrangiFilter3D(Vnew, opt);
% T = T*12./max(T(:));
% % figure, imshow3D(T);
% % keyboard;

opt.FrangiScaleRange = [1,r];
opt.FrangiScaleRatio = 1;
opt.FrangiBetaOne = .7;
opt.FrangiBetaTwo = 25;
opt.BlackWhite = false;
opt.verbose = false;
T = zeros(size(V));
scale = zeros(size(V));
for i = 1:size(V,3)
    [I,sc,~] = FrangiFilter2D(V(:,:,i), opt);
    T(:,:,i) = I;
    scale(:,:,i) = sc;
end
% figure, imshow3D(T);
T = T*12./max(T(:));
% keyboard;
end

