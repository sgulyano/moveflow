function [result, zflips, spline] = ffd_gui_both(V, swcs, track_time, soma_pos, leftpos, rightpos, dx)
%FFD_GUI_BOTH create GUI for tracking neurons both ddaD and ddaE using FFD 
%over X-axis

result = cell(size(V,3),2);         % FFD configuration
zflips = cell(size(V,3),2);         % FFD configuration
%% user parameter
if nargin < 7
    dx = 4;
end
states = 0:pi/8:pi/2-pi/8;%[0:pi/8:pi/2-pi/8 -(pi/8:pi/8:pi/2-pi/8)];%0:pi/8:pi/2-pi/8;
sel_idx = 1;
%% Init plane
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
zflip = {ones(1,ctrlpnt_dist(1)), ones(1,ctrlpnt_dist(2))};
[Ox_bs, Oz_bs] = ffd_init_from_config(config, somaX, dx, zflip);

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
f = figure('Visible','off','Position',[100,50,1050,565],'Tag','MainGUI',...
        'NumberTitle','off', 'Name','Fitting Plane to A Neuron Pair');

%  Construct the components.
uicontrol(f, 'Style','pushbutton','String','NEXT','FontSize',12,...
        'Units','normalized','Position',[0.2,0.02,0.1,0.05],...
        'Callback',@nextbutton_Callback);
uicontrol(f, 'Style','pushbutton','String','DONE','FontSize',12,...
        'Units','normalized','Position',[0.85,0.02,0.1,0.05],...
        'Callback',@donebutton_Callback);
uicontrol(f, 'Style','text','FontSize',10,...
        'Units','normalized','Position',[0.33,0.01,0.5,0.06],...
        'String',['Drag GREEN/RED points to the tips of LEFT/RIGHT dendrites '...
        'to fit the neuron. Press NEXT button (or Space) to go to the next frame.']);
ha = axes(f, 'Units','normalized','Position',[0.05,0.1,0.93,0.88]);
% draw image
himg = imshow(padarray(V(:,:,track_time+1), [20 0], 0, 'post')); title(['Time : ' num2str(track_time) '/' num2str(size(V,3)-1)]);
hold on;
% draw dendrite mask
green = cat(3, zeros(sizeI), ones(sizeI), zeros(sizeI));
hmask = imshow(green);
hmask.Visible = 'off';
set(hmask, 'AlphaData', I_mask*.2);
% draw ffd control point
px = [Ox_bs(:) Ox_bs(:)]';
py = (ones(length(Ox_bs),1)*yrange)' + ones(2,1)*Oz_bs;
hctrl = plot(px, py, 'c');
% set(hctrl, {'color'}, num2cell(parula(length(Ox_bs)), 2));
pyt = [py(1,:)'; fliplr(py(2,:))']; pxt = [px(1,:)'; fliplr(px(2,:))'];
hctrl_plane = plot(pxt([1:end 1]), pyt([1:end 1]), 'c');
pyt = [py(1,sel_idx:sel_idx+1)'; fliplr(py(2,sel_idx:sel_idx+1))'];
pxt = [px(1,sel_idx:sel_idx+1)'; fliplr(px(2,sel_idx:sel_idx+1))'];
hsel_stripe = plot(pxt([1:end 1]), pyt([1:end 1]), 'r--');
hsel_stripe.Visible = 'off';
hpntL = impoint(ha, leftpos(track_time+1,1), leftpos(track_time+1,2));
setColor(hpntL, 'r');
hpntR = impoint(ha, rightpos(track_time+1,1), rightpos(track_time+1,2));
setColor(hpntR, 'g');


hold off;
% draw point for interactive FFD fitting
N = xrange(2)-xrange(1)+1;

set(f, 'KeyPressFcn', {@guiKeyPress, N});
f.Visible = 'on';

%%
result(track_time+1,:) = config;
zflips(track_time+1,:) = zflip;

% perturb config
DONE = false;
for t = track_time+2:size(V,3)
    set(himg, 'CData', padarray(V(:,:,t), [20 0], 0, 'post')); title(['Time : ' num2str(t-1) '/' num2str(size(V,3)-1)]);
    update_trace(N)
    uiwait(gcf);
    result(t,:) = config;
    zflips(t,:) = zflip;
    if DONE, break; end
end

delete(gcf);

    

    function update_trace(N)
        % update ffd control points
        ddadX = soma_pos{1}(t,1);
        ddadY = soma_pos{1}(t,2);
        ddaeX = soma_pos{2}(t,1);
        ddaeY = soma_pos{2}(t,2);
        newsomaX = (ddadX + ddaeX)/2;
        dsomaY = round(somaY - (ddadY+ddaeY)/2);
        [Ox,Oz] = ffd_init_from_config(config, newsomaX, dx, zflip);
        px1 = [Ox(:) Ox(:)]';
        newyrange = yrange - dsomaY;
        py1 = (ones(length(Ox),1)*newyrange)' + ones(2,1)*Oz;
        set(hctrl, {'XData'}, num2cell(px1, 1)')
        set(hctrl, {'YData'}, num2cell(py1, 1)')
        pyt1 = [py1(1,:)'; fliplr(py1(2,:))']; pxt1 = [px1(1,:)'; fliplr(px1(2,:))'];
        
        
        set(hctrl_plane, 'XData', pxt1([1:end 1]), 'YData', pyt1([1:end 1]));
        pyt1 = [py1(1,sel_idx:sel_idx+1)'; fliplr(py1(2,sel_idx:sel_idx+1))'];
        pxt1 = [px1(1,sel_idx:sel_idx+1)'; fliplr(px1(2,sel_idx:sel_idx+1))'];
        set(hsel_stripe, 'XData', pxt1([1:end 1]), 'YData', pyt1([1:end 1]));
        
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
        [E,~,Es] = ffd_energy(I_bs, I, Tx, newyrange, config, ...
                [1,leftpos(t,1); ddad_somaidx, ddadX; ddae_somaidx, ddaeX; N,rightpos(t,1)]);
        fprintf('E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', E, Es{1}, Es{2}, Es{3}, Es{4});
    end

    function donebutton_Callback(~,~) 
        % Display surf plot of the currently selected data.
        DONE = true;
        uiresume(gcbf)
    end

    function guiKeyPress(~, EventData, N)
        if sel_idx > length(config{1})
            iic = 2;
            jjc = sel_idx - length(config{1});
        else
            iic = 1;
            jjc = length(config{1}) - sel_idx + 1;
        end
        
        
        switch EventData.Key
            case 'space'
                drawnow;
                uiresume(gcbf)
            case 'leftarrow'
                if sel_idx > 1
                    sel_idx = sel_idx - 1;
                end
            case 'rightarrow'
                if sel_idx < length(config{1})+length(config{2})
                    sel_idx = sel_idx + 1;
                end
            case 'w'
                zflip{iic}(jjc) = -1;
            case 's'
                zflip{iic}(jjc) = 1;
            case 'a'
                ist = find(config{iic}(jjc) == states);
                ist = ist - 1;
                if ist > 0
                    config{iic}(jjc) = states(ist);
                end
            case 'd'
                ist = find(config{iic}(jjc) == states);
                ist = ist + 1;
                if ist <= length(states)
                    config{iic}(jjc) = states(ist);
                end
        end
        update_trace(N);
    end
end

function nextbutton_Callback(~,~) 
    % Display surf plot of the currently selected data.
    drawnow;
    uiresume(gcbf)
end



