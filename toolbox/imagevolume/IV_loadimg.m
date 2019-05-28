function V = IV_loadimg( folder )
%LOADIMAGE load image into 3D array

filelist = dir([folder '*']);
cur = 1;
V = [];
for f = filelist'
    C = strsplit(f.name, '.');
    if f.name(1) == '.' || f.isdir || ~ismember(C{end}, {'tif'})
        continue
    end
    X = imread([folder f.name], C{end});
    % convert to grayscale
    if ndims(X) == 3
        X = rgb2gray(X);
    end
    X = single(X);
    % normalize image
    X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
    V(:,:,cur) = X;
    cur = cur + 1;
end

if isempty(V)
    warning(['ERROR: no image volume at ' folder])
end
end

