function [result, leftpos, rightpos, spline] = ffd_gui(V, swc, track_time, soma_pos, dx)
%FFD_GUI create GUI for tracking neuron using FFD over X-axis

result = cell(size(V,3),2);         % FFD configuration
leftpos = zeros(size(V,3),2);       % left-side pos of tips of dendrites
rightpos = zeros(size(V,3),2);      % right-side pos of tips of dendrites
%% user parameter
if nargin < 5
    dx = 4;
end
states = 0:pi/8:pi/2-pi/8;

%%
% convert swc to image
sizeI = [size(V,1), size(V,2)];
% create baseline mask
I_mask = swc2pixel( swc, sizeI );
I_mask = imdilate(I_mask, ones(3));

% crop image and relocate soma
[ii, jj] = find(I_mask);
% get bounding box
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
% locate soma
somaX = round(swc(1,3));
somaY = round(swc(1,4));
xx = xrange(1):xrange(2);
somaidx = find(xx == somaX);
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

%% draw figure
f = figure('Visible','off','Position',[160,200,1050,785],'Tag','MainGUI');
set(f, 'KeyPressFcn', @guiKeyPress);

%  Construct the components.
uicontrol(f, 'Style','pushbutton','String','Next',...
        'Units','normalized','Position',[0.87,0.4,0.1,0.05],...
        'Callback',@nextbutton_Callback);
uicontrol(f, 'Style','pushbutton','String','Done',...
        'Units','normalized','Position',[0.87,0.55,0.1,0.05],...
        'Callback',@donebutton_Callback);
ha = axes(f, 'Units','normalized','Position',[0.05,0.1,0.8,0.8]); 

himg = imshow(V(:,:,track_time+1)); title(['Time : ' num2str(track_time)]);
hold on;
% draw dendrite mask
green = cat(3, zeros(sizeI), ones(sizeI), zeros(sizeI));
hmask = imshow(green);
set(hmask, 'AlphaData', I_mask*.2);
% draw ffd control point
px = [Ox_bs(:) Ox_bs(:)]';
py = (ones(length(Ox_bs),1)*yrange)' + ones(2,1)*Oz_bs;
hctrl = plot(px, py);
set(hctrl, {'color'}, num2cell(parula(length(Ox_bs)), 2));
hold off;
% draw point for interactive FFD fitting
N = xrange(2)-xrange(1)+1;
hpntL = impoint(ha, xrange(1), somaY);
addNewPositionCallback(hpntL, @(x)pnt_Callback(x,1,N));
fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim')+[.5 -.5], get(gca,'YLim')+[.5 -.5]);
setPositionConstraintFcn(hpntL, fcn);            % Enforce boundary constraint function
setColor(hpntL, 'r');

hpntR = impoint(ha, xrange(2), somaY);
addNewPositionCallback(hpntR, @(x)pnt_Callback(x,0,N));
setPositionConstraintFcn(hpntR, fcn);            % Enforce boundary constraint function
setColor(hpntR, 'g');
f.Visible = 'on';

%%
result(track_time+1,:) = config;
leftpos(track_time+1,:) = [xrange(1), somaY];
rightpos(track_time+1,:) = [xrange(2), somaY];

% perturb config
DONE = false;
for t = track_time+2:size(V,3)
    leftpos(t,:) = leftpos(t-1,:) + [soma_pos(t,1) - soma_pos(t-1,1) 0];
    rightpos(t,:) = rightpos(t-1,:) + [soma_pos(t,1) - soma_pos(t-1,1) 0];
    set(himg, 'CData', V(:,:,t)); title(['Time : ' num2str(t-1)]);
    
    hpntL.setPosition(leftpos(t,1),leftpos(t,2))
    hpntR.setPosition(rightpos(t,1),rightpos(t,2))

    uiwait(gcf);
    result(t,:) = config;
    disp(config)
    if DONE, break; end
end

delete(gcf);

    function pnt_Callback(h,isleft,N)
        I = im2double(V(:,:,t));
        % update leftpot and Ox position
        if isleft
            leftpos(t,:) = h;
        else
            rightpos(t,:) = h;
        end
        for iter = 1:5
            config_prev = config;
            is_change = false;
            % perturb ffd
%             keyboard;
            for i = 1:length(config)
                for j = 1:length(config{i})
                    E = inf(1,length(states));
                    for k = 1:length(states)
                        config{i}(j) = states(k);
                        Ox = ffd_init_from_config(config, soma_pos(t,1), dx);
                        Tx = ffd_interpolate(Ox, spline );
                        E(k) = ffd_energy(I_bs, I, Tx, yrange - round(somaY - soma_pos(t,2)), ...
                                config, [1,leftpos(t,1); somaidx,soma_pos(t,1); N,rightpos(t,1)]);
                    end
                    [~,pos] = min(E);
                    config{i}(j) = states(pos);
                    if states(pos) ~= config_prev{i}(j)
                        is_change = true;
                    end
                end
            end
            if ~is_change, break; end
        end
        disp(h)
        update_trace()
    end

    function update_trace()
        % update ffd control points
        [Ox,Oz] = ffd_init_from_config(config, soma_pos(t,1), dx);
        px1 = [Ox(:) Ox(:)]';
        newyrange = yrange - round(somaY - soma_pos(t,2));
        py1 = (ones(length(Ox),1)*newyrange)' + ones(2,1)*Oz;
        set(hctrl, {'XData'}, num2cell(px1, 1)')
        set(hctrl, {'YData'}, num2cell(py1, 1)')
        % update dendrite mask
        Tx = ffd_interpolate(Ox, spline);
        [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
        cx = cx(:); cy = cy(:);
        idx = I_mask(yrange(1):yrange(2),xrange(1):xrange(2));
        I_mask_t = zeros(sizeI);
        idxIn = cy >= 1 & cy <= sizeI(1) & cx >= 1 & cx <= sizeI(2);
        I_mask_t(sub2ind(sizeI, cy(idx(:)&idxIn), cx(idx(:)&idxIn))) = 1;
        set(hmask, 'AlphaData', I_mask_t*.2);
        %
        I = im2double(V(:,:,t));
        [E,~,Es] = ffd_energy(I_bs, I, Tx, newyrange, ...
                config, [1,leftpos(t,1); somaidx,soma_pos(t,1)]);
        fprintf('E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', E, Es{1}, Es{2}, Es{3}, Es{4});
    end

    function donebutton_Callback(~,~) 
        % Display surf plot of the currently selected data.
        DONE = true;
        uiresume(gcbf)
    end
end

function nextbutton_Callback(~,~) 
    % Display surf plot of the currently selected data.
    drawnow;
    uiresume(gcbf)
end

function guiKeyPress(~, EventData, ~)
    if strcmp(EventData.Key, 'space')
        drawnow;
        uiresume(gcbf)
    end
end



