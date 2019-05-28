function [Vrm, Iout, scale] = preprocess_road( tube, varargin )
%PREPROCESS_ROAD preprocess step for aerial road images

%% parameters
defaultoptions = struct('sig_curv', 20, 's_csf', 4, 'gau_val', 0, 'r_oof', 4, ...
        'responsetype', 4, 'intensity_thr', 50, 'region_size', 500);
options = get_options(defaultoptions, varargin);

t = tube;
t(t < 0) = 0; t(t > 255) = 255;
Vrm = round(t);

opts.responsetype = options.responsetype;
Vpad = padarray(Vrm, options.r_oof+3*ones(1,3), 'both');
[Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);

scale = max(scale, 3);


end

