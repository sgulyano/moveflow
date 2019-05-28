function [ I ] = align_slices( V, pickslice, xrange, yrange, cen )
%ALIGN_SLICES align slices of V with respect to pickframes by considering
%the regions within xrange and yrange centered at init_cen.

if all(cen > 0)
    s = size(V);
    ran0raw = bsxfun(@plus, cen, [xrange' yrange']);
    ran0 = [floor(ran0raw(1,:)); ceil(ran0raw(2,:))];
    ran0 = min(max(ran0, 1), ones(2,1)*s([2 1]));
    % down the stack
    I0 = V(ran0(1,2):ran0(2,2), ran0(1,1):ran0(2,1), pickslice+1);
    for z = pickslice+2:size(V,3)
%         Iz = padarray(V(:,:,z), [RANGE RANGE], 0, 'both');
        ncc = normxcorr2(I0, V(:,:,z));
        
        [~, xpeak] = find(ncc==max(ncc(:)));
        xoffSet = xpeak-size(I0,2);
        V(:,:,z) = imtranslate(V(:,:,z), [ran0(1,1) - xoffSet - 1, 0]);
%         disp([ran0(1,1) - xoffSet, 0])
        
        I0 = V(ran0(1,2):ran0(2,2), ran0(1,1):ran0(2,1), z);
    end
    % up the stack
    I0 = V(ran0(1,2):ran0(2,2), ran0(1,1):ran0(2,1), pickslice+1);
    for z = pickslice:-1:1
%         Iz = padarray(V(:,:,z), [RANGE RANGE], 0, 'both');
        ncc = normxcorr2(I0, V(:,:,z));
        
        [~, xpeak] = find(ncc==max(ncc(:)));
        xoffSet = xpeak-size(I0,2);
        V(:,:,z) = imtranslate(V(:,:,z), [ran0(1,1) - xoffSet - 1, 0]);
%         disp([ran0(1,1) - xoffSet, 0]);
        
        I0 = V(ran0(1,2):ran0(2,2), ran0(1,1):ran0(2,1), z);
    end
end
I = max(V,[],3);
end

