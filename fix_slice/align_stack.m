function [ Vnew, Valign ] = align_stack( V, Vori, refid )
DEBUG = false;
SCALE = 0.5;
%ALIGN_STACK align image stack along Z-axis
Vnew = zeros(size(V));
Vnew(:,:,refid) = im2double(V(:,:,refid));
Valign = Vori;

% Type of registration error used see registration_error.m
options.type='cc';

% Use fast forward instead of central error gradient
options.centralgrad=false;

% Use cubic interpolation
options.interpolation='linear';

% b-spline grid spacing in x and y direction
Spacing = [size(V,1) 64];

% Optimizer parameters
optim=struct('Display','off','GradObj','on','HessUpdate','lbfgs','MaxIter',30,'DiffMinChange',0.01,'DiffMaxChange',1,'TolFun',1e-16,'StoreN',5,'GoalsExactAchieve',0);

% downward
for z = refid+1:size(V,3)
    I1 = im2double(imresize(V(:,:,z),SCALE));           % target
    I2 = imresize(Vnew(:,:,z-1),SCALE);                 % reference

    % Make the Initial b-spline registration grid
    O_trans = make_init_grid(Spacing*SCALE,size(I1));

    % Convert all values tot type double
    O_trans = double(O_trans); 
    
    % Reshape O_trans from a matrix to a vector.
    sizes = size(O_trans); O_trans = O_trans(:);
    
    % Start the b-spline nonrigid registration optimizer
    O_trans = fminlbfgs(@(x)bspline_registration_gradient(x,sizes,Spacing*SCALE,I1,I2,options),O_trans,optim);

    % Reshape O_trans from a vector to a matrix
    O_trans = reshape(O_trans/SCALE,sizes);

    % Transform the input image with the found optimal grid.
    Vnew(:,:,z) = bspline_transform(O_trans,im2double(V(:,:,z)),Spacing,3);
    Valign(:,:,z) = im2uint8(bspline_transform(O_trans,im2double(Vori(:,:,z)),Spacing,3));
    
    if DEBUG
        % Make a (transformed) grid image
        Igrid = make_grid_image(Spacing,size(I1));
        Igrid=bspline_transform(O_trans,Igrid,Spacing); 
        % Show the registration results
        figure(z),
        set(gcf, 'Position', [100 100 1300 1000]);
        C = imfuse(I1,I2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        subplot(3,1,1), imshow(C); title('Green-1 Red-2');
        C = imfuse(Icor,I2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        subplot(3,1,2), imshow(C); title('new img');
        subplot(3,1,3), imshow(Igrid); title('grid');
        drawnow;
    end
end

% upward
for z = refid-1:-1:1
    I1 = im2double(imresize(V(:,:,z),SCALE));           % target
    I2 = imresize(Vnew(:,:,z+1),SCALE);                 % reference
    
    % Make the Initial b-spline registration grid
    O_trans = make_init_grid(Spacing*SCALE,size(I1));

    % Convert all values tot type double
    O_trans = double(O_trans); 
    
    % Reshape O_trans from a matrix to a vector.
    sizes = size(O_trans); O_trans = O_trans(:);

    % Start the b-spline nonrigid registration optimizer
    O_trans = fminlbfgs(@(x)bspline_registration_gradient(x,sizes,Spacing*SCALE,I1,I2,options),O_trans,optim);

    % Reshape O_trans from a vector to a matrix
    O_trans = reshape(O_trans/SCALE,sizes);

    % Transform the input image with the found optimal grid.
    Vnew(:,:,z) = bspline_transform(O_trans,im2double(V(:,:,z)),Spacing,3); 
    Valign(:,:,z) = im2uint8(bspline_transform(O_trans,im2double(Vori(:,:,z)),Spacing,3));
    
    if DEBUG
        % Make a (transformed) grid image
        Igrid = make_grid_image(Spacing,size(I1));
        Igrid=bspline_transform(O_trans,Igrid,Spacing); 
        % Show the registration results
        figure(z),
        set(gcf, 'Position', [100 100 1300 1000]);
        C = imfuse(I1,I2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        subplot(3,1,1), imshow(C); title('Green-1 Red-2');
        C = imfuse(Icor,I2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        subplot(3,1,2), imshow(C); title('new img');
        subplot(3,1,3), imshow(Igrid); title('grid');
        drawnow;
    end
end

if DEBUG
    figure, imshow(max(Vnew,[],3));
    keyboard;
end
end

