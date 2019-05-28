close all; clear all; clc;

addpath(genpath('Full_Flow_Source_Code'))
addpath(genpath('toolbox'))
addpath(genpath('UGM'))

%% input
filename = '~/Desktop/LipingsData/ZstackL1_3_2grayscale/larva3_2_z%d_t%03d.tif';
savefile = 'larva3_frames_report/larva3';
savename = '~/Desktop/snakegraphmodel/swc_dan/ZstackL1_3_2grayscale/larva3_2_t%03d.swc';

%% load traces
load([savefile '_fullflow.mat'])

A = cell(1,6);
A{2} = load([savefile '_track.mat']);
A{3} = load([savefile '_track_wo_dynamic.mat']);
A{4} = load([savefile '_track_wo_transition.mat']);
A{5} = load([savefile '_track_wo_repulsive.mat']);
A{6} = load([savefile '_track_wo_shape.mat']);
str = {'GT', 'All', 'No dynamic', 'No transition', 'No repuls', 'No shape'};

% read gold standard
golddir = '~/Desktop/LipingsData/GoldStandard/ZstackL1_3_2grayscale/';  % directory containing the gold standard
for t = 224:226
    A{1}.traces{t+1,1} = read_swc_file( sprintf('%slarva3_n1_t%03d.swc', golddir, t) );
end

%% plot tracked traces
offset = 0;
for t = 224:226
    I = im2double(gfpadj(:,:,t+1));
    MAXDIST = 15;

    ii = round(A{2}.traces{t+1,1}(:,4));
    jj = round(A{2}.traces{t+1,1}(:,3));
    xmin = max(min(jj)-MAXDIST, 1); ymin = max(min(ii)-MAXDIST, 1);
    xmax = min(max(jj)+MAXDIST, size(I,2)); ymax = min(max(ii)+MAXDIST, size(I,1));
    Icrop = I(ymin:ymax, xmin:xmax);

    figure(1);
    for i = 1:length(str)
        subplot(3,length(str),i+offset); imshow(Icrop); 
        hold on; EV_plot_img( [], bsxfun(@minus, A{i}.traces{t+1,1}, [0 0 xmin-1 ymin-1 0 0 0]), 'r' );
        if t == 224,
            title(str{i});
        end
        if i == 1,
            ylabel(['Time ' num2str(t)]);
        end
    end
    offset = offset + length(str);
end

%% plot data features
branch = im2double(rgb2gray(imread('branch.png')));

Options.BlackWhite = false;
Options.FrangiScaleRange = [1 3];
Options.FrangiScaleRatio = .5;
Options.FrangiBetaOne = 0.1;
Options.verbose = false;
[Ivessel2,~,~,Ibr2]=FrangiFilter2D(branch*255,Options);

figure(2); set(gcf,'Color','w');
subplot(1,3,1); imshow(branch(30:end-30,30:end-30));
subplot(1,3,2); imagesc(Ivessel2(30:end-30,30:end-30)); colormap parula; axis off; axis equal;
subplot(1,3,3); imagesc(Ibr2(30:end-30,30:end-30)); axis off; axis equal; colorbar;


%%
gfpfile = '~/Desktop/LipingsData/20170616L7_align/20170616L7_align_z%d_t%03d.tif';     % Ca images
savefile = '20170616L7/20170616L7';
num_slice = 9;

swclist = dir([savefile '_ddad_t*.swc']);
track_time = arrayfun(@(x)regexp(x.name,'.*_t(\d+).swc','tokens'), swclist);
track_time = cellfun(@(x)str2double(x{1}), track_time);

gfpinfo = imfinfo(sprintf(gfpfile,0,0));
V = zeros(gfpinfo.Height, gfpinfo.Width, num_slice+1, 'uint8');
for z = 0:num_slice
    X = imread( sprintf(gfpfile,z,track_time) );
    X = single(X);
    X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));

    n = 3;
    Idouble = im2double(X);
    avg = mean2(Idouble);
    sigma = std2(Idouble);
    X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);

    V(:,:,z+1) = X;
end

swc_ddad = read_swc_file( sprintf('%s_ddad_t%03d.swc', savefile, track_time) );
swc_ddad(:,3:5) = swc_ddad(:,3:5)+1;
swc_ddae = read_swc_file( sprintf('%s_ddae_t%03d.swc', savefile, track_time) );
swc_ddae(:,3:5) = swc_ddae(:,3:5)+1;

sizeI = size(V);

figure; clf; set(gcf,'Color','w');
[xx, yy, zz] = meshgrid(1:sizeI(2), 1:sizeI(1), -10);
hslice = surface(xx-1,yy-1,zz,max(V,[],3),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
colormap(gray(256));
set(hslice, 'AmbientStrength', 1.0)
set(hslice, 'SpecularStrength', 0.0)
set(gca,'YDir','reverse');
zlim([-10 15])

[xx, yy, zz] = meshgrid(1:sizeI(2), 1:sizeI(1), -2:12);
soma_mask = sqrt((xx - swc_ddad(1,3)).^2 + (yy - swc_ddad(1,4)).^2 + (zz - swc_ddad(1,5)).^2) < swc_ddad(1,6);
I_mask = padarray(swc2pixel( swc_ddad, sizeI ), [0 0 7], 0, 'both');
I_mask = imdilate(I_mask, ones(3,3,5)) | soma_mask;
sf = isosurface(I_mask, .5);
hmodel = patch(sf);
hmodel.FaceColor = 'red';
hmodel.EdgeColor = 'none';
set(hmodel,'SpecularStrength', 0.5)

soma_mask = sqrt((xx - swc_ddae(1,3)).^2 + (yy - swc_ddae(1,4)).^2 + (zz - swc_ddae(1,5)).^2) < swc_ddae(1,6);
I_mask = padarray(swc2pixel( swc_ddae, sizeI ), [0 0 7], 0, 'both');
I_mask = imdilate(I_mask, ones(3,3,5)) | soma_mask;
sf = isosurface(I_mask, .5);
hmodel2 = patch(sf);
hmodel2.FaceColor = 'green';
hmodel2.EdgeColor = 'none';
daspect([1 1 1])
axis tight
axis off
% box on
camlight 
lighting gouraud

camproj('perspective');
camorbit(20,0,'data',[1 0 0])
camdolly(0,0,.9);