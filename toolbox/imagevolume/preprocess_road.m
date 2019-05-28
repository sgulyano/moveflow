function [Vrm, Iout, scale] = preprocess_road( V, varargin )
%PREPROCESS_ROAD preprocess step for aerial road images

%% parameters
defaultoptions = struct('sig_curv', 20, 's_csf', 4, 'gau_val', 3, 'r_oof', 4, ...
        'responsetype', 4, 'intensity_thr', 50, 'region_size', 500);
options = get_options(defaultoptions, varargin);

Vm = V(:,:,3);
Vadj = imadjust(Vm,[0.2 0.5]);
V_csf = CSF_cos(Vadj, 1:options.s_csf);

Vrm = double(V_csf);
Vsmall = xor(bwareaopen(Vrm>options.intensity_thr,0),  bwareaopen(Vrm>options.intensity_thr, options.region_size));
Vrm(Vsmall) = 0;

B = medfilt2(Vrm);
B = B * 255 ./ max(B(:));

enhI = myvesselEnhanceLDE(B);
enhI = enhI * 255 ./ max(enhI(:));

opts.responsetype = options.responsetype;
Vpad = padarray(enhI, options.r_oof+3*ones(1,3), 'both');
[Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);

scale = max(scale, 2);

end

