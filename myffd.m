addpath('fminlbfgs_version2c')
addpath(genpath('Full_Flow_Source_Code'))
addpath('utils')

dataset = 1;
method  = 'smooth';
switch dataset
    case 1
        savefile = 'larva3';
        savename = '~/Desktop/snakegraphmodel/swc_dan/ZstackL1_3_2grayscale/larva3_2_t%03d.swc';
    case 2
        savefile = '1122larva1';
        savename = '~/Desktop/snakegraphmodel/swc_dan/1122Larva1/1122Larva1_t%03d.swc';
    case 3
        savefile = '1202DLa2';
        savename = '~/Desktop/snakegraphmodel/swc_dan/1202DLa2/1202DLa2_t%03d.swc';
end

load(['results/' savefile '_fullflow.mat'])

Spacing = [16 16]; % must be square
Stepsize = 0.1;
SAVE_VIDEO = true;
vector_space = 8;

%%
switch method
    case 'ffd'
        Isize = size(U(:,:,1));

        ffd_x = zeros((size(U,1)/Spacing(1))+4, (size(U,2)/Spacing(2))+4, size(U,3));
        ffd_y = zeros((size(U,1)/Spacing(1))+4, (size(U,2)/Spacing(2))+4, size(U,3));
        errx = zeros(size(U,3),1); erry = zeros(size(U,3),1); err = zeros(size(U,3),1);

        % initialize cubic B-spline grid
        O_trans = init_grid(Spacing, Isize);
        Osize = size(O_trans);
        % init Bu (which equals to Bv as well)
        u = (0:Spacing(1)-1)'/Spacing(1);
        Bu = [(1-u).^3, 3*u.^3-6*u.^2+4, -3*u.^3+3*u.^2+3*u+1, u.^3]./6;
        parfor t = 1:size(U,3)
            % get optical flow
            u_flow = U(:,:,t);
            v_flow = V(:,:,t);

            % solve O_grad by BFGS optimization
            % based on dD(q;x) = sum_{k=0}^3 sum_{l=0}^3 B_k(u) B_l(v) dF_{i+k,j+l}
            options = optimset('GradObj','on');
            O_transx_init = padarray(u_flow(Spacing(1):Spacing(1):end,Spacing(2):Spacing(2):end), [2 2], 'replicate', 'both');
            O_transy_init = padarray(v_flow(Spacing(1):Spacing(1):end,Spacing(2):Spacing(2):end), [2 2], 'replicate', 'both');
            O_gradx = fminlbfgs(@(x)bspl_ffd(x, Bu, u_flow, Spacing, Stepsize), zeros(Osize(1:2)), options);
            O_grady = fminlbfgs(@(x)bspl_ffd(x, Bu, v_flow, Spacing, Stepsize), zeros(Osize(1:2)), options);

            ffd_x(:,:,t) = O_trans(:,:,2) + O_gradx;
            ffd_y(:,:,t) = O_trans(:,:,1) + O_grady;

            % get results
            u_vec_move = zeros(Isize); v_vec_move = zeros(Isize);
            i = 1;
            for y = 1:Spacing(1):Isize(1)
                j = 1;
                for x = 1:Spacing(2):Isize(2)
                    dFx = O_gradx(i:i+3, j:j+3);
                    dFy = O_grady(i:i+3, j:j+3);

                    u_vec_move(y:y+Spacing(1)-1, x:x+Spacing(2)-1) = Bu*dFx*Bu';
                    v_vec_move(y:y+Spacing(1)-1, x:x+Spacing(2)-1) = Bu*dFy*Bu';

                    j = j+1;
                end
                i = i+1;
            end

            % check error
            errx(t) = sum(sum(abs(u_vec_move - u_flow).^2));
            erry(t) = sum(sum(abs(v_vec_move - v_flow).^2));
            err(t) = sum(sum((u_vec_move - u_flow).^2 + (v_vec_move - v_flow).^2));
            fprintf('Time step %d => Err X: %f, Err Y: %f, Err: %f\n', t, errx(t), erry(t), err(t));
        end

        figure(2), plot(err); xlabel('time step'); ylabel('total vector size square error');
        save(['results/' savefile '_ffd.mat'], 'ffd_x', 'ffd_y', 'errx', 'erry', 'err', 'Spacing', '-v7.3');
    

        %% record results in video
        if SAVE_VIDEO,
            outputVideo = VideoWriter(['results/' savefile '_ffd_color.avi']);
            open(outputVideo);

            O_trans = init_grid(Spacing, Isize);
            h1 = [];
            for t = 1:size(U,3)
                figure(1); 
                % plot ffd color
                plot_ffd_color( Is(:,:,t+1), O_trans(:,:,2), O_trans(:,:,1), ffd_x(:,:,t), ffd_y(:,:,t) )

        %         % plot ffd grid
        %         plot_ffd_grid( Is(:,:,t), Is(:,:,t+1), O_trans(:,:,2), O_trans(:,:,1), ffd_x(:,:,t), ffd_y(:,:,t) )

        %         % plot original vs ffd flow
        %         u_vec_move = zeros(Isize); v_vec_move = zeros(Isize);
        %         i = 1;
        %         for y = 1:Spacing(1):Isize(1)
        %             j = 1;
        %             for x = 1:Spacing(2):Isize(2)
        %                 dFx = ffd_x(i:i+3, j:j+3,t) - O_trans(i:i+3, j:j+3, 2);
        %                 dFy = ffd_y(i:i+3, j:j+3,t) - O_trans(i:i+3, j:j+3, 1);
        % 
        %                 u_vec_move(y:y+Spacing(1)-1, x:x+Spacing(2)-1) = Bu*dFx*Bu';
        %                 v_vec_move(y:y+Spacing(1)-1, x:x+Spacing(2)-1) = Bu*dFy*Bu';
        % 
        %                 j = j+1;
        %             end
        %             i = i+1;
        %         end
        %         if isempty(h1)
        %             subplot(1,2,1); h1 = quiver(U(:,:,t), V(:,:,t)); title('FullFlow');
        %             subplot(1,2,2); h2 = quiver(u_vec_move, v_vec_move);
        %         else
        %             set(h1, 'UData', U(:,:,t));
        %             set(h1, 'VData', V(:,:,t));
        %             set(h2, 'UData', u_vec_move);
        %             set(h2, 'VData', v_vec_move);
        %         end
        %         
        %         title(['Fulflow + FFD ' num2str(t)]);
        %         drawnow;

                set(gcf, 'Position', [0 0 800 400]);
                writeVideo(outputVideo, getframe(gcf));
            end
            close(outputVideo); 
        end
        
        
    case 'smooth'
        Isize = size(U(:,:,1));

        O_trans = init_grid(Spacing, Isize);
        Osize = size(O_trans); Osize = Osize(1:2);
        sop_x = zeros([Osize, size(U,3)]);
        sop_y = zeros([Osize, size(U,3)]);

        dUdx = diff(U,[],2); dUdx(:,end+1,:) = 0;
        dUdy = diff(U,[],1); dUdy(end+1,:,:) = 0;
        
        disp('Time Step: ');
        for t = 1:size(U,3)
            fprintf('.');
            if mod(t,100) == 0
                fprintf('\n');
            end
            % get optical flow
            u_flow = dUdx(:,:,t);%U(:,:,t);
            v_flow = dUdy(:,:,t);%V(:,:,t);

            for i = 3:size(O_trans,1)-1
                for j = 3:size(O_trans,2)-1
                    if O_trans(i,j,1) <= Isize(1) && O_trans(i,j,2) <= Isize(2)
                        sop_x(i,j,t) = mean(mean(u_flow(O_trans(i-1,j,1)+1:O_trans(i,j,1), O_trans(i,j-1,2)+1:O_trans(i,j,2))));
                        sop_y(i,j,t) = mean(mean(v_flow(O_trans(i-1,j,1)+1:O_trans(i,j,1), O_trans(i,j-1,2)+1:O_trans(i,j,2))));
                    end
                end
            end
        end
        fprintf('\n');
        
        % record in video
        if SAVE_VIDEO,
            outputVideo = VideoWriter(['results/' savefile '_smooth_ffd_color.avi']);
            open(outputVideo);
        end
        
        for t = 1:size(U,3)
            figure(1); plot_ffd_color( Is(:,:,t), O_trans(:,:,2), O_trans(:,:,1), O_trans(:,:,2)+sop_x(:,:,t), O_trans(:,:,1)+sop_y(:,:,t), 0.5 );
            title(num2str(t));
            drawnow;
            if SAVE_VIDEO,
                set(gcf, 'Position', [0 0 1200 600]);
                writeVideo(outputVideo, getframe(gcf));
            end
        end
        
        if SAVE_VIDEO, close(outputVideo); end
end


%% plot vector field
% plot result and record in video
if SAVE_VIDEO,
    outputVideo = VideoWriter(['results/' savefile '_fullflow_trace.avi']);
    open(outputVideo);
    
    vecmask = false(size(U,1), size(U,2));
    vecmask(vector_space:vector_space:size(U,1),vector_space:vector_space:size(U,2)) = true;
    h1 = [];
    hs1 = [];
    for t = 1:size(U,3)
        if ~isempty(hs1),
            for h = hs1', delete(h); end;
        end
        % read trace
        swc = read_swc_file( sprintf(savename, t) );
        swc(:,3:4) = swc(:,3:4)+1;
        
        sc = max(max(sqrt(U(:,:,t).^2+V(:,:,t).^2)));
        figure(2);
        if isempty(h1) || ~isvalid(h1)
            [h1, h2] = plotFlow(U(:,:,t), V(:,:,t), Is(:,:,t), vector_space, sc/2);
        else
            set(h1, 'CData', Is(:,:,t));
            Utmp = U(:,:,t); Utmp(~vecmask) = 0;
            set(h2, 'UData', Utmp);
            Vtmp = V(:,:,t); Vtmp(~vecmask) = 0;
            set(h2, 'VData', Vtmp);
            set(h2, 'AutoScaleFactor', sc/2);
        end
        title(num2str(t));
        
        hs1 = EV_plot_img( [], swc );
        
        set(gcf, 'Position', [0 0 800 400]);
        writeVideo(outputVideo, getframe(gcf));
    end
    close(outputVideo);
end

%% plot color coded vector field with colormap according to 
% http://members.shaw.ca/quadibloc/other/colint.htm. 
% White means no movement (saturated color). 
% Yellow = Up, Green+Blue = Left, Red = Right, Blue+Purple = Down.

% rad = sqrt(U.^2+V.^2);
% maxrad = max(rad(:))+eps;

% Unorm = U ./ maxrad;
% Vnorm = V ./ maxrad;
dUdx = diff(U,[],2); dUdx(:,end+1,:) = 0;
dUdy = diff(U,[],1); dUdy(end+1,:,:) = 0;

maxrad = 5;
dUdx = max(min(dUdx, maxrad), -maxrad);
dUdy = max(min(dUdy, maxrad), -maxrad);

if SAVE_VIDEO,
    outputVideo = VideoWriter(['results/' savefile '_fullflow_colorcode.avi']);
    open(outputVideo);
    
    for t = 1:size(U,3)
        I = repmat(Is(:,:,t), [1 1 3]);
        
        img = computeColor(dUdx(:,:,t), dUdy(:,:,t));
        img = double(img) ./ 255;
        
        figure(3); imshow((img+I)./2 );
        title(num2str(t));
        drawnow;
        
        set(gcf, 'Position', [0 0 800 400]);
        pause(0.05);
        writeVideo(outputVideo, getframe(gcf));
    end
    close(outputVideo);
end


