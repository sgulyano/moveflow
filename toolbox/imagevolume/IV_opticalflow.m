function T = IV_opticalflow( V )
%IV_OPTICALFLOW Summary of this function goes here
%   Detailed explanation goes here

T = zeros(size(V));
keyboard;
opts.BlockSize   = 8;
opts.SearchLimit = 10;
s = size(T);

for i = 1:size(V,3)-1
    tic
    [MVx, MVy] = Bidirectional_ME(V(:,:,i), V(:,:,i+1), opts);
    toc

    % Motion Compensation
    tic
    imgMC = ScalarField(V(:,:,i), MVx, MVy);
    toc;
    T(:,:,i) = imresize(imgMC, s(1:2), 'bilinear');
end

end

