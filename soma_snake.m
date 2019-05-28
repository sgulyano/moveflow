function [ contour, centroid, flag ] = soma_snake( contour, U, V, I, centroid, direction )
%SOMA_SNAKE detect soma using snake and optical flow

% find optical flow vector of centroid
flag = false;
if direction > 0
    dx = interp2(U, contour(:,1), contour(:,2));
    dy = interp2(V, contour(:,1), contour(:,2));
else
    [xx, yy] = meshgrid(1:size(U,2), 1:size(U,1));
    % since U,V are forward mapping, do interpolate scattered data
    nx = xx + U;
    ny = yy + V;

    dx = griddata(nx, ny, -U, contour(:,1), contour(:,2));
    dy = griddata(nx, ny, -V, contour(:,1), contour(:,2));
end
if any(isnan(dx)) || any(isnan(dy))         % centroid out-of-bound
    flag = true;
    return;
end

N = size(contour,1);
% move every point based on centroid
contour = contour + [dx dy];%ones(N,1)*[dx dy];

s = size(U);
maxsize = ones(N,1)*s([2 1]);
idxOut = any(contour < 1 | contour > maxsize);
if any(idxOut)
    flag = true;
    return;
end

Iadj = imadjust(double(I)./255, [0 .5], [0 1]);

h = fspecial('gaussian',[3 3],3);
f = imfilter(Iadj,h);
[fx,fy] = AM_gradient(max(f,.3));
mag = fx.^2+fy.^2;
magnorm = mag ./ max(mag(:));

% keyboard;
alpha = .5;
beta = .5;
tau = .5;
SNAKE_ITER = 5;
SNAKE_ITER1 = 10;

Fext = AM_GVF(magnorm, .2, 10, 1);

% K = AM_VFK(2, 16, 'power',2);
% Fext = AM_VFC(magnorm, K, 1);

dc = bsxfun(@minus, contour, mean(contour));
dc = dc ./ (sqrt(sum(dc.^2,2))*ones(1,2));

vert = contour + dc*3;
vert0 = vert;
% disp('Deforming the snake ...')
for i=1:SNAKE_ITER1,
    vert = AC_deform_close(vert,alpha,beta,tau,Fext,SNAKE_ITER);
    vert = AC_remesh_close(vert,1);
    
%     figure(3); imshow(Iadj)
%     AC_display(vert,'close','r');
%     AC_display(vert0,'close','b--');
%     pause
end

if isempty(vert), keyboard; end
centroid = mean(vert);
contour = vert;


% figure(3); imshow(Iadj);
% hold on; plot(contour(:,1),contour(:,2),'r-o'); hold off;
% keyboard;

% figure(4); imshow(Iadj);
% hold on; quiver(Fext(:,:,1),Fext(:,:,2),'r'); hold off;



end

