function V = readImgSeq( filepattern, num_img, num_slice, opt)
%READIMGSEQ read image sequence
%   Detailed explanation goes here

if nargin < 4;  opt = struct();  end;
if ~isfield(opt,'switchTZ');        opt.switchTZ    = true;     end;
if ~isfield(opt,'original');        opt.original    = true;     end;
if ~isfield(opt,'sigmax');          opt.sigmax      = 5;        end;

[idir,ifile,iext] = fileparts(filepattern);
Vinfo = imfinfo(fullfile(idir, [sprintf(ifile, 0, 0) iext]));
V = zeros(Vinfo.Height, Vinfo.Width, num_img+1, 'uint8');

for t = 0:num_img
    Vt = zeros(Vinfo.Height, Vinfo.Width, num_slice+1, 'uint8');
    for z = 0:num_slice
        if opt.switchTZ
            X = imread( fullfile(idir, [sprintf(ifile, z, t) iext]) );
        else
            X = imread( fullfile(idir, [sprintf(ifile, t, z) iext]) );
        end
        if ~opt.original
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));

            n = opt.sigmax;
            Idouble = im2double(X);
            avg = mean2(Idouble);
            sigma = std2(Idouble);
            X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
        end
        Vt(:,:,z+1) = X;
    end
    V(:,:,t+1) = max(Vt,[],3);
    
    fprintf('.');
    if mod(t+1,100) == 0, fprintf('\n'); end;
end
fprintf('\n');

end

