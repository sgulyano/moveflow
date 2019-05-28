function [enhI] = myvesselEnhanceLDE(I)
%VESSELENHANCELDE Enhance the vessels using ridge detector along VFC based
%on the idea of LDE (local directional evidence)

vessel_type = 'bright'; %-- 'bright' or 'dark' vessel
I     = double(I);      %-- original signal
theta = 0:15:179;       %-- orientation space for detectors    
sigma = 1:3;        %-- scale space

if strcmp(vessel_type,'dark')
    sgn = 1;
else
    sgn = -1;
end
if ismatrix(I)
    [nr,nc] = size(I);
    stack = zeros(nr,nc,length(sigma));
    for ii = 1 : length(sigma)
        stack(:,:,ii) = LDE2d(I,theta,sigma(ii),sgn);
    end
    enhI = max(stack,[],3);
else
    error('image must be 2D')
end
end


function [resp] = LDE2d(I,all_theta,sigma,sgn)
%LDE Perform local directional evidence filtering
%   all_theta -- orientation for detection filter
%   psi       -- orientations of evidence filter, for each value of theta
%   d         -- d = k*sigma
drawtemplate        = 0;
l_theta             = length(all_theta);

siz     = round(4*sigma);
[X,Y]   = meshgrid(-siz:siz);
G       = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
Gxx     = (X.^2-sigma^2).*G/(sigma^4);
Gyy     = (Y.^2-sigma^2).*G/(sigma^4);
Gxy     = (X.*Y).*G/(sigma^4);



% compute filter response
resp = zeros(size(I));
for iixy = 1 : l_theta
    rotXY = all_theta(iixy);
    
    % compute filter
    v = [-cosd(rotXY), sind(rotXY)];
    v = v ./ sqrt(sum(v.^2));
    R_d = Gxx*v(1)^2 + Gyy*v(2)^2 + Gxy*v(1)*v(2);

    % compute response
    I_theta = convnfft(I,sgn*sigma^1.5*R_d,'same');

    % combind nearby response
    resp = max(resp, I_theta);

    if drawtemplate
        figure(2);
        keyboard;
        imagesc(R_d(:,:,13));
        drawnow; colormap('hot'); axis off;
        patch(isocaps(R_d+R_b+R_f,0.01),...
            'FaceColor','interp',...
            'EdgeColor','none')
    end

end
end



