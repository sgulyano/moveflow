load('myfullflow.mat')

Spacing = [16 16]; % must be square
SAVE_VIDEO = true;

Isize = size(U(:,:,1));

O_trans = init_grid(Spacing, Isize);
Osize = size(O_trans); Osize = Osize(1:2);
sop_x = zeros([Osize, size(U,3)]);
sop_y = zeros([Osize, size(U,3)]);

% record in video
if SAVE_VIDEO,
    outputVideo = VideoWriter('demo_smooth_flow.avi');
    open(outputVideo);
end

h1 = [];
for t = 1:size(U,3)
    disp(['Time Step: ', num2str(t)]);
    % get optical flow
    u_flow = U(:,:,t);
    v_flow = V(:,:,t);

    for i = 3:size(O_trans,1)-1
        for j = 3:size(O_trans,2)-1
            if O_trans(i,j,1) <= Isize(1) && O_trans(i,j,2) <= Isize(2)
                sop_x(i,j,t) = mean(mean(u_flow(O_trans(i-1,j,1)+1:O_trans(i,j,1), O_trans(i,j-1,2)+1:O_trans(i,j,2))));
                sop_y(i,j,t) = mean(mean(v_flow(O_trans(i-1,j,1)+1:O_trans(i,j,1), O_trans(i,j-1,2)+1:O_trans(i,j,2))));
            end
        end
    end
    figure(1); plot_ffd_color( Is(:,:,t), O_trans(:,:,2), O_trans(:,:,1), O_trans(:,:,2)+sop_x(:,:,t), O_trans(:,:,1)+sop_y(:,:,t) )
    title(num2str(t));
    
    if SAVE_VIDEO,
        set(gcf, 'Position', [0 0 1200 600]);
        writeVideo(outputVideo, getframe(gcf));
    end
end

if SAVE_VIDEO, close(outputVideo); end
