function [Vsmall, Vrm, Iout, scale, offset] = IV_preprocess_flycircuit( V, ZRATIO, varargin )
%IV_PREPROCESS_OP preprocess step for image volume of OP from DIADEM

%% parameters
defaultoptions = struct('s_csf', 4, 'r_oof', 4, 'responsetype', 4);
options = get_options(defaultoptions, varargin);

% crop image
[yy,xx,zz] = ind2sub(size(V), find(V>1));
offset = [min(xx), min(yy), min(zz)];
Vsmall = V(min(yy):max(yy), min(xx):max(xx), min(zz):max(zz));

FNUM = size(Vsmall,3);

figure(1); imshow(max(Vsmall,[],3),[]);
%%
disp('Preprocess for OP Dataset');
V_adj = zeros(size(Vsmall));
for fcount = 1:FNUM
    V_adj(:,:,fcount) = imadjust(uint8(Vsmall(:,:,fcount)),[0,0.3],[0,1]);
end
%figure(2); imshow3D(V_adj);

V_csf = CSF_cos(uint8(V_adj), 1:options.s_csf);
%figure(3); imshow3D(V_csf);
Vrm = double(V_csf);

%% OOF
if nargout > 2
    opts.spacing = [1, 1, ZRATIO];
    opts.responsetype = options.responsetype;
    Vpad = padarray(Vrm, options.r_oof+3*ones(1,3), 'both');
    [Iout, scale] = oof3response(Vpad, 1:options.r_oof, opts);
    figure(7), imshow3D(Iout);
    
    scale = max(scale,2);
end
end