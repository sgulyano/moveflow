function [ snks, flag ] = snaketrack( snks, U, V, I, opt )
%SNAKETRACK tracking neuron using combination of snake and optical flow
% update snake points

if nargin < 5;  opt = struct();  end;
if ~isfield(opt,'alpha');       opt.alpha     = 1.0;    end;
if ~isfield(opt,'beta');        opt.beta      = 1.0;    end;
if ~isfield(opt,'tau');         opt.tau       = 1.0;    end;
if ~isfield(opt,'oof_thr');     opt.oof_thr   = 2;      end;
if ~isfield(opt,'r_oof');       opt.r_oof     = 4;      end;
if ~isfield(opt,'sigmas');      opt.sigmas    = 2;      end;
if ~isfield(opt,'iter');        opt.iter      = 20;     end;
if ~isfield(opt,'remesh_it');   opt.remesh_it = 4;      end;
if ~isfield(opt,'OF_thr');      opt.OF_thr = 50;        end;
if ~isfield(opt,'DEBUG');       opt.DEBUG     = false;  end;

I = double(I); I = I ./ max(I(:)); I = I * 255;

count_out = 0;
for i = 1:length(snks)
    dx = interp2(U, snks(i).vert(:,1), snks(i).vert(:,2));
    dy = interp2(V, snks(i).vert(:,1), snks(i).vert(:,2));

    dx(isnan(dx)) = 0;
    dy(isnan(dy)) = 0;
    
    % init snake by optical flow
    n = size(snks(i).vert,1);
    snks(i).vert = snks(i).vert + [dx dy zeros(n,1)];
    snks(i) = AC_remesh(snks(i));

    count_out = count_out + sum(any( isnan(snks(i).vert),2 ));
end
if opt.DEBUG
    figure(14); EV_plot_img( I./255, snake2swc(snks) ); title('optical flow'); drawnow;
end

if mean(mean(U.^2+V.^2)) < opt.OF_thr
    disp(' : use only optical flow');
else
    disp(' : use snake initialized by optical flow');
    
    % preprocessing
    mu = .2;
    GVF_ITER = 10;
    h = fspecial('gaussian',[3 3],3);
    f = imfilter(I,h);
    Fext = AM_GVF(f, mu, GVF_ITER, 1);
    Fexts = {Fext(:,:,1), Fext(:,:,2)};

    opts.spacing = [1, 1, 2];
    opts.responsetype = 4;
    Vpad = padarray(I, opt.r_oof+3*ones(1,3), 'both');
    Iout = oof3response(Vpad, 1:opt.r_oof, opts);

    [Vx, Vy, ~] = IV_eig(Iout, opt.sigmas);
    ev = {Vx, Vy};

    mask = Iout>opt.oof_thr;

    s = [size(I) 1];

    if opt.DEBUG
        figure(12), imshow(I,[]);
        % hold on; quiver(ev{1}, ev{2}); hold off;
        hold on; quiver(Fexts{1}, Fexts{2}); hold off;
    end
    
    count_out = 0;
    for i = 1:length(snks)
        % evolve snake
        for itt = 1:opt.iter
            snks(i) = snake_deform(snks(i),opt.alpha,opt.beta,opt.tau,Fexts,ev,mask);
            
            if opt.DEBUG
                figure(11), imshow(I,[]);
                hold on;
                contour(mask, [.5 .5], 'g');
                for j = 1:length(snks)
                    if j == i
                        plot(snks(i).vert(:,1), snks(i).vert(:,2), 'r-', 'LineWidth', 2); 
                    else
                        plot(snks(j).vert(:,1), snks(j).vert(:,2), 'b', 'LineWidth', 2);
                    end
                end
                hold off; drawnow;
                title([i itt]);
                keyboard;
            end
            
            if mod(itt, opt.remesh_it) == 0
                snks(i) = AC_remesh(snks(i));
            end
        end

    %     % ensure branch point coherency
    %     if snks(i).collide(1) > 0
    %         snks(i).vert(1,:) = snks(snks(i).collide(1)).vert(snks(i).collide(2),:);
    %     end

        count_out = count_out + sum(any( isnan(snks(i).vert),2 ));
    end
end
% check if half of neuron is out-of-bound, then tracking is done
total_len = sum(arrayfun(@(x)size(x.vert,1), snks));
if count_out > total_len/2
    flag = true;
else
    flag = false;
end

