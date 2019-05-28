function test_main_exp_combine()
clear all
% close all
clc

rng(48)

addpath(genpath('toolbox/utils'))
addpath(genpath('toolbox/minf'))
addpath(genpath('fix_slice'))

dataset = 1;

switch dataset
    case 1
        directory = 'C:\Users\y_sar\Documents\MATLAB\DanData\1110DLa_2_P1';
        filename = '1110DLa_2_P1_z%d_t%02d.tif';
        savedir = 'C:\Users\y_sar\Documents\moveflow\1110DLa_2_P1';
        mi_range = 13:23;
    case 2
        directory = 'C:\Users\y_sar\Documents\MATLAB\DanData\GcamP images\WT\0322L1_100_P1';
        filename = '0322L1_100_P1_z%d_t%03d.tif';
        savedir = 'C:\Users\y_sar\Documents\moveflow\0322L1_100_P1';
        mi_range = 38:68;
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
    zflips = cell(size(V,3),2);         % FFD configuration
    leftpos = dd_plane(i).leftpos;       % left-side pos of tips of dendrites
    rightpos = dd_plane(i).rightpos;     % right-side pos of tips of dendrites
    % track parameter
    dx = 4;
    states = 0:pi/8:pi/2-pi/8;%[0:pi/8:pi/2-pi/8 -(pi/8:pi/8:pi/2-pi/8)];%0:pi/8:pi/2-pi/8;

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
    zflip = {ones(size(config{1})), ones(size(config{2}))};
    zflips(track_time+1,:) = zflip;

    Oz_prev = Oz_bs;
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

        Scurrent = config; Zcurrent = zflip;
        [Ox, Oz] = ffd_init_from_config(Scurrent, newsomaX, dx, Zcurrent);
        Tx = ffd_interpolate(Ox, spline );
        [Ecurrent1,~,Es] = ffd_energy(I_bs, I, Tx, yrange - dsomaY, Scurrent, ...
                [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N, rightpos(t,1)]);
        Ecurrent2 = ffd_depth_energy(Scurrent, Zcurrent, Oz_prev);
        Ecurrent = Ecurrent1 + Ecurrent2;
%         fprintf('Init : E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', Ecurrent, Es{1}, Es{2}, Es{3}, Es{4});
        Sbest = Scurrent; Ebest = Ecurrent;
        Zbest = Zcurrent; Zbest = Zcurrent;
        
        maxiter = 100;
        max_temp = 10;
        temp_change = 0.95;
        temp = max_temp;
        NctrlpntL = length(config{1});
        Nctrlpnt = NctrlpntL + length(config{2});
        for iter = 1:maxiter
            temp = temp * temp_change;
            Snew = Scurrent;
            Znew = Zcurrent;
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
                for iz = [-1 1]
                new_state_i = state_i + dstate_i;
                if new_state_i < 1 || new_state_i > length(states)
                    continue;
                end
                
                Snew{gr}(ii) = states(new_state_i);
                Znew{gr}(ii) = iz;
%                 if ii > NctrlpntL
%                     Snew{2}(ii-NctrlpntL) = states(new_state_i);
%                 else
%                     Snew{1}(ii) = states(new_state_i);
%                 end

                [Ox, Oz] = ffd_init_from_config(Snew, newsomaX, dx, Znew);
                Tx = ffd_interpolate(Ox, spline );
                Enew1 = ffd_energy(I_bs, I, Tx, yrange - dsomaY, Snew, ...
                        [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N, rightpos(t,1)]);
                Enew2 = ffd_depth_energy(Snew, Znew, Oz_prev);
                Enew = Enew1 + Enew2;
                if Enew <= Ecurrent || exp( (Ecurrent-Enew)/temp ) > rand()
                    Scurrent = Snew; Zcurrent = Znew;
                    Ecurrent = Enew;
                end
                if Enew <= Ebest
                    Sbest = Snew; Zbest = Znew;
                    Ebest = Enew;
                end
                end
            end
        end
        config = Sbest;
        zflip = Zbest;

        % update ffd control points
        [Ox,Oz] = ffd_init_from_config(config, newsomaX, dx, zflip);
        newyrange = yrange - dsomaY;

        % update dendrite mask
        Tx = ffd_interpolate(Ox, spline);
        I = im2double(V(:,:,t));
        [E,~,Es] = ffd_energy(I_bs, I, Tx, newyrange, config, ...
                [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N,rightpos(t,1)]);
%         fprintf('E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', E, Es{1}, Es{2}, Es{3}, Es{4});
        configs(t,:) = config;
        zflips(t,:) = zflip;
        time_elapse(t) = toc;
    end
    
    dd_plane(i).configs = configs;
    dd_plane(i).spline = spline;
end

opt.directory = directory;
opt.filename = filename;
opt.zflip = zflips;
for i = 1:num_pairs
    mi_score = plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)], opt);
end
fprintf('Mean MI : %f', mean(mi_score(mi_range)));
fprintf(',   STD MI : %f\n', std(mi_score(mi_range)));
disp('Done');
keyboard;



function [E, Es] = ffd_depth_energy(config, S, Oz_prev)
ang1 = cumsum(config{1});
ang2 = cumsum(config{2});

Ox = {-cos(ang1)*dx, cos(ang2)*dx};
Oz = {abs(sin(ang1)*dx), abs(sin(ang2)*dx)};
    
Ox_all = [fliplr(cumsum(Ox{1})) 0 cumsum(Ox{2})];
Oz_all = [fliplr(cumsum( Oz{1}.*S{1}  )) 0 cumsum( Oz{2}.*S{2} )];

% Edep = sum(exp(Oz_all - min(Oz_all)));
Edep = sum((Oz_all - min(Oz_all)).^2);
Ecur = sum(diff(Ox_all,2).^2 + diff(Oz_all,2).^2);
Eprev = sum((Oz_all - Oz_prev).^2);

lines = [Ox_all(1:end-1)' Oz_all(1:end-1)' Ox_all(2:end)' Oz_all(2:end)'];
intrsect = false(size(lines,1));
% for i = size(lines,1)-11:size(lines,1)
for it = 1:size(lines,1)
    intrsect(it,:) = line_seg_intersect(lines(it,:), lines);
end

if any(any(tril(intrsect,-2) | triu(intrsect,2)))
    % self-intersection occurs
    E = inf;
    Es = [inf inf inf];
else
    E = 0*Edep/2 + Ecur/4 + Eprev/4;
    Es = [0*Edep/2 Ecur/4 Eprev/4];
end
end

function out = line_seg_intersect(l1, l2)
% https://en.wikipedia.org/wiki/Intersection_(Euclidean_geometry)
tol = 0.5;
float_tol = 1e-10;

dx1 = l1(:,3) - l1(:,1);
dy1 = l1(:,4) - l1(:,2);
dx2 = l2(:,3) - l2(:,1);
dy2 = l2(:,4) - l2(:,2);
cx = l2(:,1) - l1(:,1);
cy = l2(:,2) - l1(:,2);

parallel_idx = abs(dx1*dy2 - dx2*dy1) < float_tol;
s = (cy.*dx2 - cx.*dy2) ./ (dx2.*dy1 - dx1.*dy2);
tt = (dx1*cy - dy1*cx) ./ (dx2.*dy1 - dx1.*dy2);
% for non-parallel
out1 = isfinite(s) & isfinite(tt) & s >= 0 & s <= 1 & tt >= 0 & tt <= 1;
% for parallel
p1 = (cx*dx1 + cy*dy1) / sqrt(dx1^2+dy1^2);
d1 = sqrt(cx.^2 + cy.^2 - p1.^2);
p1abs = p1 / sqrt(dx1^2+dy1^2);

cx2 = l2(:,3) - l1(:,1);
cy2 = l2(:,4) - l1(:,2);
p2 = (cx2*dx1 + cy2*dy1) / sqrt(dx1^2+dy1^2);
d2 = sqrt(cx2.^2 + cy2.^2 - p2.^2);
p2abs = p2 / sqrt(dx1^2+dy1^2);

out2 = parallel_idx & (d1 <= tol | d2 <= tol) & ((p1abs >= -float_tol & p1abs <= 1+float_tol) | (p2abs >= -float_tol & p2abs <= 1+float_tol));
out = out1 | out2;
end
end