function O_trans = init_grid(Spacing, Isize)
% Make the Initial b-spline registration grid
if length(Spacing) == 3
     % Determine grid spacing
    dx=Spacing(1); dy=Spacing(2); dz=Spacing(3);
    
    % Calculate te grid coordinates (make the grid)
    [X,Y,Z]=ndgrid(-dx:dx:(Isize(1)+(dx*2)),-dy:dy:(Isize(2)+(dy*2)),-dz:dz:(Isize(3)+(dz*2)));
    O_trans=ones(size(X,1),size(X,2),size(X,3),3);
    O_trans(:,:,:,1)=X;
    O_trans(:,:,:,2)=Y;
    O_trans(:,:,:,3)=Z;
elseif length(Spacing) == 2
    % Determine grid spacing
    dx=Spacing(1); dy=Spacing(2);

    % Calculate te grid coordinates (make the grid)
    [X,Y]=ndgrid(-dx:dx:(Isize(1)+(dx*2)),-dy:dy:(Isize(2)+(dy*2)));
    O_trans=ones(size(X,1),size(X,2),2);
    O_trans(:,:,1)=X;
    O_trans(:,:,2)=Y;
else
    error('Unsupported datatype');
end
