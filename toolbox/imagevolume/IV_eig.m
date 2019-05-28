function [Vx, Vy, Vz] = IV_eig( I, sigmas )
%IV_EIG find eigenvector along neurite

% Use single or double for calculations
if(~isa(I,'double')), I=single(I); end

if ismatrix(I)
    % Calculate 2D hessian
    [Dxx, Dxy, Dyy] = Hessian2D(I,sigmas);

    if(sigmas>0)
        % Correct for scaling
        c=(sigmas^2);
        Dxx = c*Dxx;
        Dxy = c*Dxy;
        Dyy = c*Dyy;
    end

    % Calculate eigen values
    [~,~,Vx,Vy]=eig2image(Dxx,Dxy,Dyy);
    Vz = zeros(size(Vx));    
%     [~,~,Direction] = FrangiFilter2D(I);
%     Vx = sin(Direction);
%     Vy = cos(Direction);
%     Vz = zeros(size(Direction));
else
    % Calculate 3D hessian
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigmas);

    if(sigmas>0)
        % Correct for scaling
        c=(sigmas^2);
        Dxx = c*Dxx; Dxy = c*Dxy;
        Dxz = c*Dxz; Dyy = c*Dyy;
        Dyz = c*Dyz; Dzz = c*Dzz;
    end

    % Calculate eigen values
    [~,~,~,Vy,Vx,Vz]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
end

end

