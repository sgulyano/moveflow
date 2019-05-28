function energy = get_GVF_energy( f, MU, Fext )
% GET_GVF_ENERGY get GVF energy
%     Inputs
%     f           edge map
%     MU       	  a parameter for weighting edge map and mag of vec field
%     Fext        GVF field
%     
%     Outputs
%     energy      the energy map of GVF field
%

f = single(f);
fmin  = min(f(:));
fmax  = max(f(:));
f = (f-fmin)/(fmax-fmin);           % Normalize f to the range [0,1]

if ismatrix(f),
    [fx,fy] = AM_gradient(f);       % Calculate the gradient of the edge map
    fz = 0;
    u = Fext(:,:,1); v = Fext(:,:,2);
    
    [ux,uy] = AM_gradient(u);
    term1 = MU*((ux.^2) + (uy.^2));
    [vx,vy] = AM_gradient(v);
    term1 = term1 + MU*((vx.^2) + (vy.^2));
    
    tmp = (u-fx).*(u-fx) + (v-fy).*(v-fy);
else
    [fx,fy,fz] = AM_gradient(f);    % Calculate the gradient of the edge map
    u = Fext(:,:,:,1); v = Fext(:,:,:,2); w = Fext(:,:,:,3);
    
    [ux,uy,uz] = AM_gradient(u);
    term1 = MU*((ux.^2) + (uy.^2) + (uz.^2));
    [vx,vy,vz] = AM_gradient(v);
    term1 = term1 + MU*((vx.^2) + (vy.^2) + (vz.^2));
    [wx,wy,wz] = AM_gradient(w);
    term1 = term1 + MU*((wx.^2) + (wy.^2) + (wz.^2));
    
    tmp = (u-fx).*(u-fx) + (v-fy).*(v-fy) + (w-fz).*(w-fz);
end
SqrMagf = fx.*fx + fy.*fy + fz.*fz;

term2 = (SqrMagf.*tmp);

energy = (term2 + term1);

end

