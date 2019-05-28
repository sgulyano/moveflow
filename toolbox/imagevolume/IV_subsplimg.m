function [ imgdir ] = IV_subsplimg( folder, factor )
%IV_SUBSPLIMG subsampling image
imgdir = 'CFtmp/';
mkdir(imgdir);

filelist = dir([folder '*']);
cur = 1;
for f = filelist'
    C = strsplit(f.name, '.');
    if f.name(1) == '.' || f.isdir || ~ismember(C{end}, {'tif'})
        continue
    end
    im = imread([folder f.name]);
    X = fun(im, factor);
    saveimg = [imgdir, sprintf('%03d.tif', cur)];
    imwrite(X, saveimg)
%     fun = @(block_struct) {
%         imresize(block_struct.data, 1/factor)
%        };
%     blockproc([folder f.name], [400 400], @(block)fun(block, factor), 'Destination', saveimg);
    
    cur = cur + 1;
end

end


function X = fun(im, factor)
    X = imresize(im, 1/factor);
    X = double(X);
    tmp = X(:,:,2) + 2*(X(:,:,3) - X(:,:,1));
    grayIm = min(tmp,255*ones(size(tmp)));
    grayIm(grayIm<0) = 0;
    X = 255 - grayIm;
    X = X - 100;
    X(X<0) = 0;
%     keyboard;
    if any(X(:) > 255)
        keyboard;
    end
%     X = round(255*X./max(X(:)));
    X = uint8(X);
end
