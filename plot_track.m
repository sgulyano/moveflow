function mi_score = plot_track( V, dd_plane, swc_ddad, swc_ddae, soma_cen, str_title, opt )
%PLOT_TRACK plot tracking results by fitting neuron's plane

if nargin < 7;  opt = struct();  end
if ~isfield(opt,'is_gui');      opt.is_gui    = false;              end
if ~isfield(opt,'dx');          opt.dx        = 4;                  end
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end
if ~isfield(opt,'maxdep');      opt.maxdep    = 27;                 end
if ~isfield(opt,'zflip');       opt.zflip     = [];                 end
if ~isfield(opt,'directory');   opt.directory = [];                 end
if ~isfield(opt,'filename');    opt.filename  = [];                 end

if ~isempty(opt.directory) && ~isempty(opt.filename)
    [~, num_slice, switchTZ] = getNumImgAndSlice(opt.directory, opt.filename);
end

%%
col1 = [80  200 80 ]/255; % ddaD color
col2 = [235 90  235]/255; % ddaE color
% col1 = [235 90  235]/255; % ddaD color
% col2 = [80  200 80 ]/255; % ddaE color

% get starting and ending frame indices
idx = ~cellfun(@isempty, dd_plane.configs(:,1));
st_fr = find(idx, 1, 'first');
en_fr = find(idx, 1, 'last');


mi_score = zeros(1,en_fr);
dx = opt.dx;
padsize = [2 2 7];

%%
sizeI = [size(V,1), size(V,2)];
% create baseline mask
I_mask1 = swc2pixel( swc_ddad, sizeI, '2D' );
I_mask2 = swc2pixel( swc_ddae, sizeI, '2D' );
I2D_mask = imdilate(I_mask1 | I_mask2, ones(3));

% get bounding box
[ii, jj] = find(I2D_mask);
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
I_mask1_crop = I_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
I_mask2_crop = I_mask2(yrange(1):yrange(2),xrange(1):xrange(2));

ddad_somaX = round(swc_ddad(1,3));
ddad_somaY = round(swc_ddad(1,4));
ddae_somaX = round(swc_ddae(1,3));
ddae_somaY = round(swc_ddae(1,4));
somaY = round((swc_ddad(1,4) + swc_ddae(1,4))/2);

sf1 = cell(1,en_fr);
sf2 = cell(1,en_fr);
fwait = waitbar(0, ['Initialize ' str_title '... 0% complete']);
Oz_prev = [];
Esz = zeros(en_fr, 4);

xoffset1 = length(dd_plane.configs{st_fr,1})*dx;
xoffset2 = length(dd_plane.configs{st_fr,2})*dx;

for fr = st_fr:en_fr
    ddadX = round(soma_cen{1}(fr,1));
    ddadY = round(soma_cen{1}(fr,2));
    ddaeX = round(soma_cen{2}(fr,1));
    ddaeY = round(soma_cen{2}(fr,2));
    newsomaX = (ddadX + ddaeX)/2;
    newsomaY = (ddadY + ddaeY)/2;
    dsomaY = round(somaY - newsomaY);

    % init FFD
    if isempty(opt.zflip)
        [Ox,Oz] = ffd_init_from_config(dd_plane.configs(fr,:), newsomaX, dx);
        if isempty(Oz_prev)
            Oz_prev = Oz;
        end
        if sum(Oz) > 0
            [~, zflip, Esz(fr,1), Esz(fr,2:end)] = estimate_plane_z(dd_plane.configs(fr,:), dx, Oz_prev);
            [ Ox2, Oz2 ] = ffd_init_from_config( dd_plane.configs(fr,:), newsomaX, dx, zflip );
            figure(20), subplot(2,2,1); plot(Ox,Oz,'r',Ox2,Oz2,'b--'); title('Neuron Plane XZ-axis');
            subplot(2,2,2); plot(Esz); legend({'Sum', 'Dep', 'Cur', 'Pre'}); title('Break Down Z-Optim Energy');
            Oz_prev = Oz2;
            Oz = Oz2;
        end
    else
        if ismatrix(opt.zflip)
            [Ox,Oz] = ffd_init_from_config(dd_plane.configs(fr,:), newsomaX, dx, opt.zflip(fr,:));
        else
            [Ox,Oz] = ffd_init_from_config(dd_plane.configs(fr,:), newsomaX, dx, opt.zflip{fr});
        end
    end
    
    % transform dendrite in XY-plane
    Tx = ffd_interpolate(Ox, dd_plane.spline);
    newyrange = yrange - dsomaY;
    [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
    cx = cx(:); cy = cy(:);
    idxIn = cy >= 1 & cy <= sizeI(1) & cx >= 1 & cx <= sizeI(2);

    % transform dendrite in Z-axis
    Tz = ffd_interpolate(Oz, dd_plane.spline)+1;
    [cz, ~] = meshgrid(round(Tz) - min(round(Tz))+1, newyrange(1):newyrange(2));
    ddadZ = cz(ddad_somaY - yrange(1) + 1, ddad_somaX - xrange(1) + 1);
    ddaeZ = cz(ddae_somaY - yrange(1) + 1, ddae_somaX - xrange(1) + 1);
    cz = cz(:);
    sizeI3D = [sizeI 1+range(cz)];
    [xx3, yy3, zz3] = meshgrid(1:sizeI3D(2)+padsize(2)*2, 1:sizeI3D(1)+padsize(1)*2, 1:sizeI3D(3)+padsize(3)*2);
    % update 3D dendrite ddaD
    I3D = zeros(sizeI3D);
    % draw dednrite
    I3D(sub2ind(sizeI3D, cy(I_mask1_crop(:)&idxIn), cx(I_mask1_crop(:)&idxIn), cz(I_mask1_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, padsize, 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    % draw soma
    I3D((xx3 - ddadX - padsize(2)).^2 + (yy3 - ddadY - padsize(1)).^2 + 4*(zz3 - ddadZ - padsize(3)).^2 < (swc_ddad(1,6) + 2)^2) = 1;
    % get 3D isosurface
    sf1{fr} = isosurface(I3D, .5);
    sf1{fr}.vertices = bsxfun(@minus, sf1{fr}.vertices, padsize) * opt.pixelsize;
    sf1{fr}.vertices(:,3) = bsxfun(@minus, sf1{fr}.vertices(:,3), 7);
    
    mask1 = I3D(padsize(1)+1:padsize(1)+sizeI(1), padsize(2)+1:padsize(2)+sizeI(2), :);
    
    % update 3D dendrite ddaE
    I3D = zeros(sizeI3D);
    % draw dednrite
    I3D(sub2ind(sizeI3D, cy(I_mask2_crop(:)&idxIn), cx(I_mask2_crop(:)&idxIn), cz(I_mask2_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    % draw soma
    I3D((xx3 - ddaeX - padsize(2)).^2 + (yy3 - ddaeY - padsize(1)).^2 + 4*(zz3 - ddaeZ - padsize(3)).^2 < (swc_ddae(1,6) + 2)^2) = 1;    
    sf2{fr} = isosurface(I3D, .5);
    sf2{fr}.vertices = bsxfun(@minus, sf2{fr}.vertices, padsize) * opt.pixelsize;
    sf2{fr}.vertices(:,3) = bsxfun(@minus, sf2{fr}.vertices(:,3), 7);
    
    mask2 = I3D(padsize(1)+1:padsize(1)+sizeI(1), padsize(2)+1:padsize(2)+sizeI(2), :);
    
    if ~isempty(opt.directory) && ~isempty(opt.filename)
        img = zeros(sizeI(1), sizeI(2), num_slice, 'uint8');
        for z = 0:num_slice
            if switchTZ
                X = imread( fullfile(opt.directory, sprintf(opt.filename,z,fr-1)) );
            else
                X = imread( fullfile(opt.directory, sprintf(opt.filename,fr-1,z)) );
            end
            
            X = single(X);
            X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));

            n = 3;
            Idouble = im2double(X);
            avg = mean2(Idouble);
            sigma = std2(Idouble);
            X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
            
            img(:,:,z+1) = X;
        end
        masks = mask1 | mask2;
        
        % compare
        offset = ceil((padsize(3)*2 - num_slice) / 2);
        [ii,jj,kk] = meshgrid(1:sizeI(2), 1:sizeI(1), offset:offset+num_slice);
        maskint = interp3(double(masks), ii, jj, kk);
        % crop
        [yy,~,~] = ind2sub(size(img), find(maskint));
        maskcrop = imgaussfilt3(maskint(min(yy):max(yy),newsomaX-xoffset1:newsomaX+xoffset2,:), 1);
        imgcrop = img(min(yy):max(yy),newsomaX-xoffset1:newsomaX+xoffset2,:);
        figure(20);
        subplot(2,2,3); imagesc(max(maskcrop,[],3));
        subplot(2,2,4); imagesc(max(imgcrop,[],3));
%         % move img
%         zoff = padsize(3)+1 - ceil(num_slice/2);
%         if zoff > 0
%             imgmove = imtranslate(img,[0,0,zoff],'OutputView','full');
%             if size(masks,3) > size(imgmove,3)
%                 imgmove = padarray(imgmove, [0 0 size(masks,3)-size(imgmove,3)], 'post');
%             else
%                 error('wrong size');
%             end
%         else
%             error('wrong size');
%         end
%         hehe = normxcorr3(maskcrop,imgcrop, 'valid');
%         mi_score(fr) = MI2(maskcrop,imgcrop,'Normalized');
        mi_score(fr) = max(mi(maskcrop,imgcrop,64),0);
%         mi_score(fr) = mi(masks,imgmove);
        disp(['MI score : ' num2str(mi_score(fr))]);
    end
    
    waitbar(fr/en_fr, fwait, ['Initialize ' str_title '... ' num2str(100*fr/en_fr, 3) '% complete']);
end
close(fwait);

%% Construct the interface
f = figure('Visible','on','Position',[100,50,1080,580], ...
        'Name',str_title,'NumberTitle','off');

uicontrol('Style','text',...
        'Units','normalized','Position',[0.65 0.06 0.1 0.05],...
        'String','Alpha: ','FontSize',14);
txtalpha = uicontrol('Style','text','HorizontalAlignment','left',...
        'Units','normalized','Position',[0.75 0.06 0.1 0.05],...
        'String','0.5','FontSize',14);
uicontrol('Style','text',...
        'Units','normalized','Position',[0.1 0.06 0.05 0.05],...
        'String','0','FontSize',14);
uicontrol('Style','text',...
        'Units','normalized','Position',[0.55 0.06 0.05 0.05],...
        'String','1','FontSize',14);
hslider_alpha = uicontrol('Style','slider','Enable','on',...
        'Min', 0, 'Max', 1, 'Value', 0.5, ...
        'Units','normalized','Position',[0.15 0.07 0.4 0.03]);
hslider_alpha.Callback = {@alphaSlider, hslider_alpha, txtalpha};

uicontrol('Style','text',...
        'Units','normalized','Position',[0.65 0.02 0.1 0.05],...
        'String','Frame: ','FontSize',14);
txtframe = uicontrol('Style','text','HorizontalAlignment','left',...
        'Units','normalized','Position',[0.75 0.02 0.1 0.05],...
        'String',num2str(st_fr),'FontSize',14);
uicontrol('Style','text',...
        'Units','normalized','Position',[0.1 0.02 0.05 0.05],...
        'String',num2str(st_fr),'FontSize',14);
uicontrol('Style','text',...
        'Units','normalized','Position',[0.55 0.02 0.05 0.05],...
        'String',num2str(en_fr),'FontSize',14);
hslider_frame = uicontrol('Style','slider','Enable','on',...
        'Min', st_fr, 'Max', en_fr, 'Value', st_fr, 'SliderStep', [1/(en_fr-st_fr), 10/(en_fr-st_fr)],...
        'Units','normalized','Position',[0.15 0.03 0.4 0.03]);
hslider_frame.Callback = {@frameSlider, hslider_frame, txtframe};

uicontrol('Style','pushbutton',...
        'Units','normalized','Position',[0.85 0.02 0.1 0.1],...
        'String','Done','FontSize',14,...
        'Callback','uiresume(gcbf)');


ha = axes('Units','normalized','Position',[0.1,0.2,0.85,0.70],'FontSize',16);

[xx, yy, zz] = meshgrid((0:size(V,2)-1)*opt.pixelsize, (0:size(V,1)-1)*opt.pixelsize, 0);
axes(ha);
xlabel('X (\mum)');
ylabel('Y (\mum)');
zlabel('Z (\mum)');
himg = surface(xx,yy,zz,255-V(:,:,st_fr),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
colormap(gray(256));
himg.SpecularStrength = 0;
himg.DiffuseStrength = 0;
himg.AmbientStrength = 1;
alpha(himg, 0.5);

hold on;
hmodel1 = patch(sf1{st_fr});
hmodel1.FaceColor = col1;
hmodel1.EdgeColor = 'none';
% hmodel1.Visible = 'off';
hmodel2 = patch(sf2{st_fr});
hmodel2.FaceColor = col2;
hmodel2.EdgeColor = 'none';
% hmodel2.Visible = 'off';
hold off;

axis([0 sizeI(2) 0 sizeI(1)] * opt.pixelsize);
 
lighting gouraud
% view(0, 90)
view(-180, 10)
% zlim([-8 8]);
zticks([-8 -4 0 4 8])
% xlim([100 200]);
set(ha,'Ydir','reverse');
set(ha,'Zdir','reverse');
title('Time : 0 s');
box on
daspect([1 1 0.33])
camproj('perspective')
camlight;
% uiwait(f);
% close(f);

    function frameSlider(~,~,hslider,htxt)
        hslider.Value = round(hslider.Value);
        t = hslider.Value;
        set(himg, 'CData', 255-V(:,:,t));
        htxt.String = num2str(t);
        
        title(['Time : ' num2str(t*opt.tpf,3) ' s']);
        
        hmodel1.Vertices = sf1{t}.vertices;
        hmodel1.Faces = sf1{t}.faces;
        
        hmodel2.Vertices = sf2{t}.vertices;
        hmodel2.Faces = sf2{t}.faces;
    end

    function alphaSlider(~,~,hslider,htxt)
        htxt.String = num2str(hslider.Value,2);
        alpha(himg, hslider.Value);
    end

end

