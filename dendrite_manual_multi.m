function dds = dendrite_manual_multi( V, t_init, dd_init, somas )
%SOMA_MANUAL manually change soma position
t = t_init;
n_neuron = size(dd_init,1);
dds = zeros([size(dd_init) size(V,3)]);
dds(:,:,t) = dd_init;
selected = 0;
hselect = [];

% % Type of registration error used see registration_error.m
options.type='sd';
% % Use fast forward instead of central error gradient
options.centralgrad=false;
% % Use cubic interpolation
options.interpolation='cubic';
% % b-spline grid spacing in x and y direction
Spacing = [size(V,1) 32];
% 
% % Optimizer parameters
optim=struct('Display','off','GradObj','on','HessUpdate','lbfgs','MaxIter', 100,'DiffMinChange',0.01,'DiffMaxChange',1,'TolFun',1e-16,'StoreN',5,'GoalsExactAchieve',0);

sizeI = [size(V,1) size(V,2)];

fcn = makeConstrainToRectFcn('imrect',[1 size(V,2)],[1 size(V,1)]);

f = figure('Position',[10 100 1080 720], 'NumberTitle','off', ...
        'MenuBar', 'none', 'Name','Tracking Soma');
    
uimenu(f,'Label','Help','Callback',['helpdlg({''Right click to select dendrite bounding box'',',...
        '''Left click to adjust dendrite bouding box'',',...
        '''space/right arrow = next frame'',''left arrow = previous frame'',',...
        '''Press Esc when finish''})']);
ha = axes('Units','normalized','Position',[0.05,0.05,0.93,0.93],'FontSize',16);
himg = imshow(V(:,:,t), []);

hold on
hdd = cell(1,n_neuron);
for i = 1:n_neuron
    hdd{i} = plot(dd_init(i,[1 3 3 1 1]), dd_init(i, [2 2 4 4 2]), 'g--', 'LineWidth', 2);
end
hold off

set(himg,'ButtonDownFcn',@mouseCallback);

title(['Tracking Soma at Frame ' num2str(t-1)]);

set(f, 'WindowKeyPressFcn', @keyPressFcn);
uiwait(f);
pause(1);
try
    uiwait(msgbox('Done'));
catch
    warning('User hold spacebar. msgbox was closed before uiwait is activated.');
end
delete(f);


    function mouseCallback(hObject, ~)
        if strcmp('alt', f.SelectionType)
            deselect = true;
            axesHandle  = get(hObject,'Parent');
            coordinates = get(axesHandle,'CurrentPoint');
            coordinates = coordinates(1,1:2);

            for ii = 1:n_neuron
                if all(dds(ii,1:2,t) <= coordinates & coordinates <= dds(ii,3:4,t))
                    if selected > 0
                        deselectRect();
                    end
                    selected = ii;
                    hdd{selected}.XData = [];
                    hdd{selected}.YData = [];
                    hselect = imrect(ha, dds(ii,:,t) - [0 0 dds(ii,1:2,t)], 'PositionConstraintFcn', fcn);
                    deselect = false;
                    break
                end
            end
            if deselect && selected > 0
                deselectRect();
            end
        end
    end

    function keyPressFcn(~, eventdata)
        switch eventdata.Key
            case 'r'
                if ismember('control', eventdata.Modifier) && selected > 0
                    dds(selected,:,t) = 0;
                    delete(hselect);
                    selected = 0;
                end
            case {'space','rightarrow'}
                if selected > 0
                    deselectRect();
                end
                
                if t >= size(V,3)
                    uiresume;
                else
                    t = t + 1;
                    % 
                    Iprev = imgaussian(im2double(V(:,:,t-1)),1);    % moving
                    I = imgaussian(im2double(V(:,:,t)),1);          % reference

                    O_trans = make_init_grid(Spacing,sizeI);

                    % Convert all values tot type double
                    O_trans = double(O_trans); 

                    % Reshape O_trans from a matrix to a vector.
                    sizes = size(O_trans); O_trans = O_trans(:);

                    % Start the b-spline nonrigid registration optimizer
                    O_trans = fminlbfgs(@(x)bspline_registration_gradient(x,sizes,Spacing,I,Iprev,options),O_trans,optim);

                    % Reshape O_trans from a vector to a matrix
                    O_trans = reshape(O_trans,sizes);

                    % Transform the input image with the found optimal grid.
                    [~, B] = bspline_transform(O_trans,I,Spacing,3);

                    mBx = mean(B(:,:,2),1);
                    mBy = mean(B(:,:,1),2);

                    dds(:,1:2:end,t) = dds(:,1:2:end,t-1) + interp1(mBx, dds(:,1:2:end,t-1), 'linear', 0);
                    dds(:,2:2:end,t) = dds(:,2:2:end,t-1) + interp1(mBy, dds(:,2:2:end,t-1), 'linear', 0);

                    dds(:,:,t) = min(max(dds(:,:,t), 1), ones(n_neuron,1)*sizeI([2 1 2 1]));
                end
            case 'leftarrow'
                if t > 1
                    t = t - 1;
                end
            case 'escape'
                uiresume(gcbf);
        end
        title(['Tracking Soma at Frame ' num2str(t-1)]);
        if t >= 1 && t <= size(V,3),
            himg.CData = V(:,:,t);
            for ii = 1:n_neuron
                hdd{ii}.XData = dds(ii,[1 3 3 1 1],t);
                hdd{ii}.YData = dds(ii, [2 2 4 4 2],t);
            end
            pause(0.05);
        end
    end

    function deselectRect()
        pos = getPosition(hselect);
        delete(hselect);
        dds(selected,:,t) = pos + [0 0 pos(1:2)];
        hdd{selected}.XData = dds(selected,[1 3 3 1 1],t);
        hdd{selected}.YData = dds(selected, [2 2 4 4 2],t);
        selected = 0;
    end
    
end

