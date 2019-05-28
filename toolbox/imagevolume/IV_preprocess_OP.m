function [Vrm, Iout, scale] = IV_preprocess_OP( V, ZRATIO, varargin )
%IV_PREPROCESS_OP preprocess step for image volume of OP from DIADEM

%% parameters
defaultoptions = struct('sig_curv', 20, 's_csf', 2, 'gau_val', 10, 'r_oof', 8, ...
        'responsetype', 4, 'intensity_thr', 50, 'region_size', 200);
options = get_options(defaultoptions, varargin);

FNUM = size(V,3);

figure(1); imshow3D(V);
%%
disp('Preprocess for OP Dataset');
V_adj = zeros(size(V));
for fcount = 1:FNUM
    V_adj(:,:,fcount) = imadjust(uint8(V(:,:,fcount)),[0,0.3],[0,1]);
end
figure(2); imshow3D(V_adj);
Vrm = V_adj;

if options.gau_val > 0
    %% Curvelet Preprocessing
    curvelets=zeros(size(V));
    for i=1:FNUM
        [curvelets(:,:,i)]=NfdctDemo(V_adj(:,:,i), options.sig_curv);
    end


    %% LoG
    V_csf = CSF_cos(uint8(curvelets), options.s_csf);
    figure(3); imshow3D(V_csf);
    Vrm = double(V_csf);
    % %%
%     gau_val = options.gau_val;
%     B = padarray(V_csf,[gau_val, gau_val, gau_val], 'replicate', 'both');
    
%     h = fspecial('gaussian', [gau_val, gau_val], gau_val);
%     % h = fspecial('gaussian', [50, 50], 50);
%     f = imfilter(double(B),h);
%     F = f(gau_val+1:end-gau_val, gau_val+1:end-gau_val, gau_val+1:end-gau_val);
%     figure(4), imshow3D(F)
% 
%     Vsub = double(V_csf)-F;
%     Vsub(Vsub < 0) = 0;
%     Vsub = Vsub ./ max(Vsub(:)); Vsub = Vsub .* 255;
%     % figure(5), imshow3D(Vsub);
% 
%     % Iout = xor(bwareaopen(Vsub>50,0),  bwareaopen(Vsub>50,400));
%     Iout = xor(bwareaopen(Vsub>options.intensity_thr,0),  bwareaopen(Vsub>options.intensity_thr, options.region_size));
%     Vrm = Vsub;
%     Vrm(Iout) = 0;
%     figure(6), imshow3D(Vrm);
end

%% OOF
if nargout > 1
    opts.spacing = [1, 1, ZRATIO];
    opts.responsetype = options.responsetype;
    Vpad = padarray(Vrm, options.r_oof+3*ones(1,3), 'both');
    [Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);
    figure(7), imshow3D(Iout);
end
% keyboard;
end