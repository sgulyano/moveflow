function [Vrm, Iout, scale, mask] = preprocess_drive(V, mask, varargin)
% PREPROCESS_DRIVE preprocess step for retina images in DRIVE

% %% parameters
defaultoptions = struct('r_oof', 4, 'responsetype', 4);
options = get_options(defaultoptions, varargin);
% 
% % keyboard;
% 
% Vin = 255-rgb2gray(V);
% Vc = imadjust(Vin,[0.3 0.7]);
% 
% % figure, imshow(Vc)
% % %% Curvelet Preprocessing
% % curvelets=zeros(size(V));
% % for i=1:FNUM
% %     [curvelets(:,:,i)]=NfdctDemo(V(:,:,i), options.sig_curv);
% % end
% % if REPORT
% %     tmp = imcrop(mean(curvelets,3), rect);
% %     figure(2), imshow(tmp, [], 'Border','tight');
% % end
% % 
% % %% LoG
% % V_csf = zeros(size(Vin));
% % for s_csf = 1:4
% %     v_csf = CSF_cos(Vin, s_csf);
% %     V_csf = max(double(v_csf), V_csf);
% % end
% % figure, imshow(v_csf,[])
% V_csf = CSF_cos(Vc, 1:options.s_csf);
% 
% 
% % 
% % %%
% % gau_val = options.gau_val;
% % B = padarray(V_csf,[gau_val, gau_val, gau_val], 'replicate', 'both');
% % 
% % h = fspecial('gaussian', [gau_val, gau_val], gau_val);
% % % h = fspecial('gaussian', [50, 50], 50);
% % f = imfilter(double(B),h);
% % F = f(gau_val+1:end-gau_val, gau_val+1:end-gau_val, gau_val+1:end-gau_val);
% % % figure(4), imshow3D(F)
% % 
% % Vsub = double(V_csf)-F;
% % Vsub(Vsub < 0) = 0;
% % Vsub = Vsub ./ max(Vsub(:)); Vsub = Vsub .* 255;
% % % figure(5), imshow3D(Vsub);
% % 
% % % Iout = xor(bwareaopen(Vsub>50,0),  bwareaopen(Vsub>50,400));
% % Iout = xor(bwareaopen(Vsub>options.intensity_thr,0),  bwareaopen(Vsub>options.intensity_thr, options.region_size));
% % Vrm = Vsub;
% % Vrm(Iout) = 0;
% % % figure(6), imshow3D(Vrm);
% % 
% % %% OOF
% 
% Vrm = double(V_csf);
% 
% opts.responsetype = options.responsetype;
% Vpad = padarray(Vrm, options.r_oof+3*ones(1,3), 'both');
% [Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);
% scale = max(scale, 2);
% 
% Iout = min(Iout * 5, 20);
% se = strel('disk',12);
% Iout(~imerode(mask, se)) = 0;
% figure(7), imshow(Iout, []);


image = double(V) ./ 255;
output = retina_extract( image );
mask = output.segmented;

se = strel('disk',5);
Vmask = imdilate(mask, se);

Vrm = output.respimage;
Vrm(~Vmask) = 0;

opts.responsetype = options.responsetype;
Vpad = padarray(Vrm, options.r_oof+3*ones(1,3), 'both');
[Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);
scale = max(scale, 2);

Iout = min(Iout * 3, 20);
figure(7), imshow(Iout, []);


end