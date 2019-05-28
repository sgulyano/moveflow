function [gfp, gfp1sl] = readImg(directory, filename, opt)
%READIMG read image stack sequence in directory with filepattern
if nargin < 3;  opt = struct();  end
if ~isfield(opt,'normalized');          opt.normalized = true;          end
if ~isfield(opt,'adj_slicewise');       opt.adj_slicewise = true;       end

% Read image given the file pattern
[num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);
pickslice = floor(num_slice/2);

if num_img <= 0 || num_slice <= 0
    gfp = [];
    return;
end

% read image
disp('Reading GFP and RFP images');
gfpinfo = imfinfo(fullfile(directory, sprintf(filename,0,0)));
gfp = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');
gfp1sl = zeros(gfpinfo.Height, gfpinfo.Width, num_img+1, 'uint8');

fwait = waitbar(0,'Loading image sequence');
for t = 0:num_img
    waitbar(t/num_img,fwait,['Loading frame ' num2str(t) '/' num2str(num_img)]);
    V = zeros(gfpinfo.Height, gfpinfo.Width, opt.sto, 'uint8');
    for z = opt.sfr-1:opt.sto-1
        if switchTZ
            X = imread( fullfile(directory, sprintf(filename,z,t)) );
        else
            X = imread( fullfile(directory, sprintf(filename,t,z)) );
        end
        if opt.normalized
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        end
        if z == pickslice
            gfp1sl(:,:,t+1) = X;
        end
        if opt.adj_slicewise
            n = 3;
            Idouble = im2double(X);
            avg = mean2(Idouble);
            sigma = std2(Idouble);
            X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
        end
        V(:,:,z+1) = X;
    end
    gfp(:,:,t+1) = max(V(:,:,opt.sfr:opt.sto),[],3);
end
close(fwait); 
end

