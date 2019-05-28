function [T, scale] = IV_tensor3D( V )
%IV_TENSOR3D Summary of this function goes here
%   Detailed explanation goes here

r = 8;

T = zeros(size(V));
scale = zeros(size(V));

for i = 1:size(V,3)
    [s,o]=encode(V(:,:,i)); %To encode the tensorized gradient of an image
%     figure
%     subplot(3,4,1); imshow(V(:,:,i),[])
%     subplot(3,4,2); imshow(s,[])
    I = T(:,:,i);
    sc = scale(:,:,i);
    for j = 1:r
        [saliency,ballness,orientation]=vote(s,o,j); %Tensor Voting with sigma=5
%         subplot(3,4,j+2); imshow(saliency,[])
        
        [a,b] = max(cat(3, I, saliency./j), [], 3);
        I = a;
        sc(b==2) = j;
    end
    T(:,:,i) = I;
    scale(:,:,i) = sc;
%     subplot(3,4,11); imshow(I,[])
%     subplot(3,4,12); imshow(sc,[])
end
T = T*12./max(T(:));

