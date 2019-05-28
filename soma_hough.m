function [ centroid, rad, flag ] = soma_hough( U, V, centers, radii, centroid )
%SOMA_HOUGH detect soma using ellipse hough transform. See soma as ellipse

% find optical flow vector of centroid
flag = false;
rad = 0;

dx = interp2(U, centroid(1), centroid(2));
dy = interp2(V, centroid(1), centroid(2));

if any(isnan(dx)) || any(isnan(dy))         % centroid out-of-bound
    flag = true;
    return;
end

centroid = centroid + [dx dy];

[ds,pos] = min(sqrt(sum(bsxfun(@minus, centers, centroid).^2,2)));
if ds > 40,
    flag = true;
    return;
end

% figure(2); subplot(2,1,2); imshow(f);
% viscircles(centers, radii,'EdgeColor','b');
% centroid_old = centroid;
centroid = centers(pos,:);
rad = radii(pos);

% figure(11); imshow(f);
% hold on; plot(centroid_old(1),centroid_old(2),'rx'); hold off;
% viscircles(centers, radii,'EdgeColor','b');
% viscircles(centroid, rad,'EdgeColor','r');
% keyboard;

% E = edge(f,'canny');
% keyboard;
% 
% % override some default parameters
% params.minMajorAxis = 9;
% params.maxMajorAxis = 10;
% 
% % note that the edge (or gradient) image is used
% bestFits = ellipseDetection(E, params);

% figure;
% image(I);
% %ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
% ellipse(bestFits(:,3),bestFits(:,4),bestFits(:,5)*pi/180,bestFits(:,1),bestFits(:,2),'k');
% 
% fprintf('Output %d best fits.\n', size(bestFits,1));
% 
% figure;
% image(I);
% %ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
% ellipse(bestFits(:,3),bestFits(:,4),bestFits(:,5)*pi/180,bestFits(:,1),bestFits(:,2),'k');


end

