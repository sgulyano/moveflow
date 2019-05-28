function [ centroid, BW, flag ] = soma_levelset( U, V, I, centroid, BW, direction )
%SOMA_LEVELSET segment soma using level set
    flag = false;
    if direction > 0
        dx = interp2(U, centroid(1), centroid(2));
        dy = interp2(V, centroid(1), centroid(2));
    else
        dx = interp2(-U, centroid(1), centroid(2));
        dy = interp2(-V, centroid(1), centroid(2));
%         [xx, yy] = meshgrid(1:size(U,2), 1:size(U,1));
%         % since U,V are forward mapping, do interpolate scattered data
%         nx = xx + U;
%         ny = yy + V;
%         
%         dx = griddata(nx, ny, -U, centroid(1), centroid(2));
%         dy = griddata(nx, ny, -V, centroid(1), centroid(2));
    end
    if isnan(dx) || isnan(dy)               % centroid out-of-bound
        flag = true;
        return;
    end
    % centroid = centroid + [dx dy];

    % init level-set by optical flow
    [ii, jj] = ind2sub(size(BW), find(BW));
    nx = round(jj + dx); ny = round(ii + dy);
    idxOut = ny < 1 | ny > size(U,1) | nx < 1 | nx > size(U,2);

    idx = sub2ind(size(BW), ny(~idxOut), nx(~idxOut));
    BW = false(size(BW));
    BW(idx) = true;
    if isempty(idx)                         % moved object out-of-bound
        flag = true;
        return;
    end
    
    BW_old = BW;
    % applya level set
    phi = double(- bwdist(BW) + bwdist(~BW) - BW + 0.5);
%     I = double(I);
    Iadj = imadjust(double(I)./255, [0 .5], [0 1]) * 255;
    mu = 100;
    g = ac_gradient_map(Iadj, 5); 
    propagation_weight = .1; GAC_weight = 1; 
    delta_t = .5; n_iters = 10; show_result = 0;
    phi = ac_hybrid_model(Iadj-mu, phi, propagation_weight, GAC_weight, g, ...
            delta_t, n_iters, show_result);
%     phi = ac_ChanVese_model(Iadj, phi, 1, 0.0001, ...
%             delta_t, n_iters, show_result);

    if isempty(phi)                         % moved object out-of-bound
        BW = false(size(BW));
        flag = true;
        return;
    else
        BW = phi >= 0;
    end
    [ii, jj] = ind2sub(size(BW), find(BW));
    centroid = [mean(jj), mean(ii)];
    
    
    
    figure(2); imshow( I, [] );
    hold on;
    plot(centroid(1), centroid(2), 'rx');
    contour(phi, [0 0], 'r', 'LineWidth', 2);
    contour(BW_old, [.5 .5], 'b', 'LineWidth', 2);
    hold off;
    drawnow;
    keyboard;
    if sum(idxOut) > length(idxOut)/2       % half object out-of-bound
        flag = true;
    end
end

