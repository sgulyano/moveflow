function [Vrm, Iout, scale] = preproc_dan(V, varargin)
% PREPROC_DAN preprocess step for 4D image from Dan

%% parameters
defaultoptions = struct('imadjustparam', [0.02 0.25],  'r_oof', 8, ...
        'responsetype', 4);
options = get_options(defaultoptions, varargin);

FNUM = size(V,3);
%% adjust intensity
V_adj = zeros(size(V));
for fcount = 1:FNUM
    V_adj(:,:,fcount) = imadjust(V(:,:,fcount),options.imadjustparam,[0,1]);
end

%% OOF
if nargout > 1
    opts.spacing = [1, 1, 4];
    opts.responsetype = options.responsetype;
    Vpad = padarray(double(V_adj), options.r_oof+3*ones(1,3), 'both');
    [Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);
    scale = max(scale, 2);
%     figure(8), imshow3D(Iout);
end
Vrm = double(V_adj);
end