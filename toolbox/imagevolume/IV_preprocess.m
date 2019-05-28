function [Vrm, Iout, scale] = IV_preprocess(V, ZRATIO, varargin)
% PREPROCESS preprocess step for image volume of drosophila dataset

%% parameters
defaultoptions = struct('sig_curv', 20, 's_csf', 2, 'gau_val', 10, 'r_oof', 8, ...
        'responsetype', 4, 'intensity_thr', 30, 'region_size', 100);
options = get_options(defaultoptions, varargin);

FNUM = size(V,3);

REPORT = false;
if REPORT
    close all
    rect = [30, 30, 200, 200];
    tmp = imcrop(mean(V,3), rect);
    figure(1), imshow(tmp, [], 'Border','tight');
end

%% Curvelet Preprocessing
curvelets=zeros(size(V));
for i=1:FNUM
    [curvelets(:,:,i)]=NfdctDemo(V(:,:,i), options.sig_curv);
end
if REPORT
    tmp = imcrop(mean(curvelets,3), rect);
    figure(2), imshow(tmp, [], 'Border','tight');
end

%% LoG
V_csf = CSF_cos(uint8(curvelets), options.s_csf);
if REPORT
    tmp = imcrop(mean(V_csf,3), rect);
    figure(3), imshow(tmp, [], 'Border','tight');
end

%%
gau_val = options.gau_val;
B = padarray(V_csf,[gau_val, gau_val, gau_val], 'replicate', 'both');

h = fspecial('gaussian', [gau_val, gau_val], gau_val);
% h = fspecial('gaussian', [50, 50], 50);
f = imfilter(double(B),h);
F = f(gau_val+1:end-gau_val, gau_val+1:end-gau_val, gau_val+1:end-gau_val);
% figure(4), imshow3D(F)

Vsub = double(V_csf)-F;
Vsub(Vsub < 0) = 0;
Vsub = Vsub ./ max(Vsub(:)); Vsub = Vsub .* 255;
% figure(5), imshow3D(Vsub);

% Iout = xor(bwareaopen(Vsub>50,0),  bwareaopen(Vsub>50,400));
Iout = xor(bwareaopen(Vsub>options.intensity_thr,0),  bwareaopen(Vsub>options.intensity_thr, options.region_size));
Vrm = Vsub;
Vrm(Iout) = 0;
% figure(6), imshow3D(Vrm);
if REPORT
    tmp = imcrop(mean(Vsub,3), rect);
    figure(4), imshow(tmp, [], 'Border','tight');
    tmp = imcrop(mean(Vrm,3), rect);
    figure(5), imshow(tmp, [], 'Border','tight');
end

%% OOF
if nargout > 1
    opts.spacing = [1, 1, ZRATIO];
    opts.responsetype = options.responsetype;
    Vpad = padarray(Vrm, options.r_oof+3*ones(1,3), 'both');
    [Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);
    figure(7), imshow3D(Iout);
end

if REPORT
    tmp = imcrop(mean(Iout,3), rect);
    figure(6), imshow(tmp, [], 'Border','tight');
%     keyboard;
end
end