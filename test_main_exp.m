clear all
% close all
clc

rng(48)

addpath(genpath('toolbox/utils'))
addpath(genpath('toolbox/minf'))
addpath(genpath('fix_slice'))

dataset = 2;%1;%

switch dataset
    case 1
        directory = 'C:\Users\sguly\Documents\MATLAB\DanData\GcamP images\WT\1110DLa_2_P1';
        filename = '1110DLa_2_P1_z%d_t%02d.tif';
        savedir = 'C:\Users\sguly\Documents\moveflow\1110DLa_2_P1';
        mi_range = 13:23;
    case 2
        directory = 'C:\Users\sguly\Documents\MATLAB\DanData\GcamP images\WT\0322L1_100_P1';
        filename = '0322L1_100_P1_z%d_t%03d.tif';
        savedir = 'C:\Users\sguly\Documents\moveflow\0322L1_100_P1';
        mi_range = 38:68;
    case 3
        directory = 'C:\Users\sguly\Documents\MATLAB\DanData\20180404mcdGFP_Curvature\20180404-L6-B';
        filename = '20180404-L6-B_z%02d_t%02d.tif';
        savedir = 'C:\Users\sguly\Documents\moveflow\20180404-L6-B';
        mi_range = 13:23;
end

load(fullfile(savedir, 'user_input.mat'), ...
        'params', 'soma_cen', 'soma_rad', 'forwardmovement', ...
        'initswctime', 'ddad_swc', 'ddae_swc', 'dd_plane');

num_pairs = length(ddad_swc);
            
% read the middle frame for tuning parameter
[num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);

[gfp, gfp1sl] = readImg(directory, filename, params);
           

%% track a pair of ddaD and ddaE by fitting a plane through folding
for i = 1:num_pairs
    dd_plane(i).configs = [];
end

V = gfp;
time_elapse = zeros(1,size(V,3));

for i = 1:num_pairs
    swcs = {ddad_swc{i}, ddae_swc{i}};
    soma_pos = {soma_cen{1,i}, soma_cen{2,i}};
    track_time = initswctime(i);
    configs = cell(size(V,3),2);         % FFD configuration
    leftpos = dd_plane(i).leftpos;       % left-side pos of tips of dendrites
    rightpos = dd_plane(i).rightpos;     % right-side pos of tips of dendrites
    % track parameter
    dx = 4;
%     states = -pi/2+pi/8:pi/8:pi/2-pi/8;%0:pi/8:pi/2-pi/8;%[0:pi/8:pi/2-pi/8 -(pi/8:pi/8:pi/2-pi/8)];%0:pi/8:pi/2-pi/8;
    states = 0:pi/8:pi/2-pi/8;%[0:pi/8:pi/2-pi/8 -(pi/8:pi/8:pi/2-pi/8)];%0:pi/8:pi/2-pi/8;
%     states = -pi/2+pi/8:pi/8:pi/2-pi/8;

    % Init plane
    % convert swc to image
    sizeI = [size(V,1), size(V,2)];
    % create baseline mask
    I_mask = swc2pixel( swcs{1}, sizeI ) | swc2pixel( swcs{2}, sizeI );
    I_mask = imdilate(I_mask, ones(3));

    % crop image and relocate soma
    [ii, jj] = find(I_mask);
    % get bounding box
    xrange = [min(jj) max(jj)];
    yrange = [min(ii) max(ii)];
    N = xrange(2)-xrange(1)+1;
    % locate soma
    ddad_somaX = round(swcs{1}(1,3));
    ddad_somaY = round(swcs{1}(1,4));
    ddae_somaX = round(swcs{2}(1,3));
    ddae_somaY = round(swcs{2}(1,4));
    somaX = round((ddad_somaX + ddae_somaX)/2);
    somaY = round((ddad_somaY + ddae_somaY)/2);
    xx = xrange(1):xrange(2);
    ddad_somaidx = find(xx == ddad_somaX);
    ddae_somaidx = find(xx == ddae_somaX);
    I_bs = im2double(V(yrange(1):yrange(2),xrange(1):xrange(2),track_time+1));

    % init FFD
    ctrlpnt_dist = ceil(abs(xrange-somaX)/dx)+[1 2];
    config = {zeros(1,ctrlpnt_dist(1)), zeros(1,ctrlpnt_dist(2))};
    [Ox_bs, Oz_bs] = ffd_init_from_config(config, somaX, dx);

    % init spline
    spline.Bu = zeros(4,dx);
    uu = 0:dx-1;
    u = (uu/dx) - floor(uu/dx);
    spline.Bu(1,:) = (1-u).^3/6;
    spline.Bu(2,:) = ( 3*u.^3 - 6*u.^2 + 4)/6;
    spline.Bu(3,:) = (-3*u.^3 + 3*u.^2 + 3*u + 1)/6;
    spline.Bu(4,:) = u.^3/6;
    spline.u_index_array = mod(xx-Ox_bs(1),dx)*4;   % times 4, it is column in Bu
    spline.i_array = floor((xx-Ox_bs(1))/dx); 

    %
    configs(track_time+1,:) = config;

    % perturb config
    for t = track_time+2:size(V,3)
        I = im2double(V(:,:,t));

        tic;
        % Optimize using Simulated Annealing
        ddadX = soma_pos{1}(t,1);
        ddadY = soma_pos{1}(t,2);
        ddaeX = soma_pos{2}(t,1);
        ddaeY = soma_pos{2}(t,2);
        newsomaX = (ddadX + ddaeX)/2;
        dsomaY = round(somaY - (ddadY+ddaeY)/2);

        Scurrent = config;
        Ox = ffd_init_from_config(Scurrent, newsomaX, dx);
        Tx = ffd_interpolate(Ox, spline );
        [Ecurrent,~,Es] = ffd_energy(I_bs, I, Tx, yrange - dsomaY, Scurrent, ...
                [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N, rightpos(t,1)]);
%         fprintf('Init : E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', Ecurrent, Es{1}, Es{2}, Es{3}, Es{4});
        Sbest = Scurrent; Ebest = Ecurrent;
        maxiter = 100;
        max_temp = 10;
        temp_change = 0.95;
        temp = max_temp;
        NctrlpntL = length(config{1});
        Nctrlpnt = NctrlpntL + length(config{2});
        for iter = 1:maxiter
            temp = temp * temp_change;
            Snew = Scurrent;
            % perturb
            gr = randi(2);
            ii = randi(length(config{gr}));
            [~, state_i] = max(states == Snew{gr}(ii));
%             ii = randi(Nctrlpnt);
%             if ii > NctrlpntL
%                 [~, state_i] = max(states == Snew{2}(ii-NctrlpntL));
%             else
%                 [~, state_i] = max(states == Snew{1}(ii));
%             end

            for dstate_i = [-1 1]
                new_state_i = state_i + dstate_i;
                if new_state_i < 1 || new_state_i > length(states)
                    continue;
                end
                
                Snew{gr}(ii) = states(new_state_i);
%                 if ii > NctrlpntL
%                     Snew{2}(ii-NctrlpntL) = states(new_state_i);
%                 else
%                     Snew{1}(ii) = states(new_state_i);
%                 end

                Ox = ffd_init_from_config(Snew, newsomaX, dx);
                Tx = ffd_interpolate(Ox, spline );
                Enew = ffd_energy(I_bs, I, Tx, yrange - dsomaY, Snew, ...
                        [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N, rightpos(t,1)]);
                if Enew <= Ecurrent || exp( (Ecurrent-Enew)/temp ) > rand()
                    Scurrent = Snew;
                    Ecurrent = Enew;
                end
                if Enew <= Ebest
                    Sbest = Snew;
                    Ebest = Enew;
                end
            end
        end
        config = Sbest;
        
        configs(t,:) = config;
        time_elapse(t) = toc;
        
        % update ffd control points
        [Ox,Oz] = ffd_init_from_config(config, newsomaX, dx);
        newyrange = yrange - dsomaY;

        % update dendrite mask
        Tx = ffd_interpolate(Ox, spline);
        I = im2double(V(:,:,t));
        [E,~,Es] = ffd_energy(I_bs, I, Tx, newyrange, config, ...
                [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N,rightpos(t,1)]);
%         fprintf('E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', E, Es{1}, Es{2}, Es{3}, Es{4});
    end
    
    dd_plane(i).configs = configs;
    dd_plane(i).spline = spline;
end
fprintf('Time Elapsed : %f s\n', mean(time_elapse));

opt.directory = directory;
opt.filename = filename;
for i = 1:num_pairs
    mi_score = plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)], opt);
end
%%
ffd_plot_exp(fullfile(savedir, 'accordion.mp4'), gfp, {ddad_swc{i}, ddae_swc{i}}, ...
        initswctime(i), {soma_cen{1,i}, soma_cen{2,i}}, ...
        dd_plane(i));
%%
fprintf('Mean MI : %f', mean(mi_score(mi_range)));
fprintf(',   STD MI : %f\n', std(mi_score(mi_range)));
fprintf('Time Elapsed : %f s\n', mean(time_elapse));
disp('Done');