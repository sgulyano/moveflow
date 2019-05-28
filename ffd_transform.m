function I_trans = ffd_transform( I_bs_crop, Tx, yrange, sizeI )
%FFD_TRANSFORM Given baseline image around neuron I_bs_crop and transformed
%coordinate Tx and new yrange of the current image with image size sizeI,
%return the transformed image

[xx, yy] = meshgrid(Tx, yrange(1):yrange(2));
[xq, yq] = meshgrid(1:sizeI(2), 1:sizeI(1));

I_trans = griddata(xx(:),yy(:),double(I_bs_crop(:)),xq,yq);
I_trans(isnan(I_trans)) = 0;

end

