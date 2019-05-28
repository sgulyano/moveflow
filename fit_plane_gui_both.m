function [result, leftpos, rightpos] = fit_plane_gui_both(V, swcs, track_time, soma_pos, dx)
%FFD_GUI_BOTH create GUI for tracking neurons both ddaD and ddaE using FFD 
%over X-axis

%% user parameter
if nargin < 5
    dx = 8;
end
states = [-1 0 1];

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
ddadX_init = round(swcs{1}(1,3));
ddadY_init = round(swcs{1}(1,4));
ddaeX_init = round(swcs{2}(1,3));
ddaeY_init = round(swcs{2}(1,4));
cenX = round((ddadX_init + ddaeX_init)/2);
cenY = round((ddadY_init + ddaeY_init)/2);
I_bs = im2double(V(yrange(1):yrange(2),xrange(1):xrange(2),track_time+1));

% init FFD
num_ctrlpnt = ceil(abs(xrange-cenX)/dx) + [1 1];
Sinit = {dx*ones(1,num_ctrlpnt(1)), 8*ones(1,num_ctrlpnt(2))};
[Ox_bs, Oz_bs] = plane_init_from_config(Sinit, cenX, dx);


%% draw figure
f = figure('Visible','off','Position',[160,50,1050,550],'Tag','MainGUI',...
        'NumberTitle','off', 'Name','Fitting Plane to A Neuron Pair');
set(f, 'KeyPressFcn', @guiKeyPress);

%  Construct the components.
uicontrol(f, 'Style','pushbutton','String','NEXT','FontSize',12,...
        'Units','normalized','Position',[0.2,0.02,0.1,0.05],...
        'Callback',@nextbutton_Callback);
uicontrol(f, 'Style','pushbutton','String','DONE','FontSize',12,...
        'Units','normalized','Position',[0.85,0.02,0.1,0.05],...
        'Callback',@donebutton_Callback);
uicontrol(f, 'Style','text','FontSize',12,...
        'Units','normalized','Position',[0.33,0.01,0.5,0.08],...
        'String',['Drag GREEN/RED points to the tips of LEFT/RIGHT dendrites '...
        'to fit the neuron. Press NEXT button (or Space) to go to the next frame.']);
ha = axes(f, 'Units','normalized','Position',[0.05,0.1,0.93,0.88]);
% draw image
himg = imshow(V(:,:,track_time+1)); title(['Time : ' num2str(track_time) '/' num2str(size(V,3)-1)]);
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
pyt = [py(1,:)'; fliplr(py(2,:))']; pxt = [px(1,:)'; fliplr(px(2,:))'];
hctrl_plane = plot(pxt([1:end 1]), pyt([1:end 1]), 'b--');
hold off;
% draw point for interactive FFD fitting
hpntL = impoint(ha, xrange(1), cenY);
addNewPositionCallback(hpntL, @(x)pnt_Callback(x,1));
fcn = makeConstrainToRectFcn('impoint', get(gca,'XLim')+[.5 -.5], get(gca,'YLim')+[.5 -.5]);
setPositionConstraintFcn(hpntL, fcn);            % Enforce boundary constraint function
setColor(hpntL, 'r');

hpntR = impoint(ha, xrange(2), cenY);
addNewPositionCallback(hpntR, @(x)pnt_Callback(x,0));
setPositionConstraintFcn(hpntR, fcn);            % Enforce boundary constraint function
setColor(hpntR, 'g');
f.Visible = 'on';

%%
result = cell(size(V,3),2);         % FFD configuration
leftpos = zeros(size(V,3),2);       % left-side pos of tips of dendrites
rightpos = zeros(size(V,3),2);      % right-side pos of tips of dendrites

result(track_time+1,:) = Sinit;
leftpos(track_time+1,:) = [xrange(1), cenY];
rightpos(track_time+1,:) = [xrange(2), cenY];

% perturb config
DONE = false;
for t = track_time+2:size(V,3)
    Sbest = Sinit;
    ds = mean(cellfun(@(x)(x(t,1) - x(t-1,1)), soma_pos));
    leftpos(t,:) = leftpos(t-1,:) + [ds 0];
    rightpos(t,:) = rightpos(t-1,:) + [ds 0];
    set(himg, 'CData', V(:,:,t)); title(['Time : ' num2str(t-1) '/' num2str(size(V,3)-1)]);

    hpntL.setPosition(leftpos(t,1),leftpos(t,2))
    hpntR.setPosition(rightpos(t,1),rightpos(t,2))
    
    uiwait(gcf);
    result(t,:) = Sbest;
    if DONE, break; end
end

delete(gcf);

    function pnt_Callback(h,isleft)
        I = im2double(V(:,:,t));
        % update leftpot and Ox position
        if isleft
            leftpos(t,:) = h;
            str = 'Left';
        else
            rightpos(t,:) = h;
            str = 'Right';
        end
        % Optimize using Simulated Annealing
        ddadX = soma_pos{1}(t,1);
        ddadY = soma_pos{1}(t,2);
        ddaeX = soma_pos{2}(t,1);
        ddaeY = soma_pos{2}(t,2);
        cenXt = (ddadX + ddaeX)/2;
        dcenY = round(cenY - (ddadY+ddaeY)/2);
        
        Scurrent = Sinit;
        Ox = plane_init_from_config(Scurrent, cenXt);
        [Ecurrent,~,Es] = fit_plane_energy(I_bs, I, Ox_bs, Ox, xrange, yrange - dcenY, dx, ...
                [xrange(1),leftpos(t,1); ddadX_init, ddadX; ddaeX_init, ddaeX; xrange(2), rightpos(t,1)]);
        fprintf('Init : E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', Ecurrent, Es{1}, Es{2}, Es{3}, Es{4});
        Sbest = Scurrent; Ebest = Ecurrent;
        maxiter = 100;
        max_temp = 100;
        temp_change = 0.95;
        temp = max_temp;
        for iter = 1:maxiter
            temp = temp * temp_change;
            Snew = Scurrent;
            % perturb
            gr = randi(2);
            i = randi(length(Sinit{gr}));
            
            for state_i = states
                pos_i = Scurrent{gr}(i) + state_i;
                if abs(Snew{gr}(i)) > dx
                    continue;
                end
                Snew{gr}(i) = pos_i;

                Ox = plane_init_from_config(Snew, cenXt);
                [Enew,~,~] = fit_plane_energy(I_bs, I, Ox_bs, Ox, xrange, yrange - dcenY, dx, ...
                        [xrange(1),leftpos(t,1); ddadX_init, ddadX; ddaeX_init, ddaeX; xrange(2), rightpos(t,1)]);
                %fprintf('New (%d) : E=%f Eimg=%f Esmooth=%f Eshape=%f Econ=%f\n', state_i, Enew, Es{1}, Es{2}, Es{3}, Es{4});
                
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

        fprintf('%s Point Position = [%.2f, %.2f]\n', str, h);
        update_trace();
    end

    function update_trace()
        % update ffd control points
        ddadX = soma_pos{1}(t,1);
        ddadY = soma_pos{1}(t,2);
        ddaeX = soma_pos{2}(t,1);
        ddaeY = soma_pos{2}(t,2);
        cenXt = (ddadX + ddaeX)/2;
        dcenY = round(cenY - (ddadY+ddaeY)/2);
        
        [Ox, Oz] = plane_init_from_config(Sbest, cenXt, dx);
        px1 = [Ox(:) Ox(:)]';
        newyrange = yrange - dcenY;
        py1 = (ones(length(Ox),1)*newyrange)' + ones(2,1)*Oz;
        set(hctrl, {'XData'}, num2cell(px1, 1)')
        set(hctrl, {'YData'}, num2cell(py1, 1)')
        pyt1 = [py1(1,:)'; fliplr(py1(2,:))']; pxt1 = [px1(1,:)'; fliplr(px1(2,:))'];
        set(hctrl_plane, 'XData', pxt1([1:end 1]), 'YData', pyt1([1:end 1]));
        
        % update dendrite mask
        Tx = interp1(Ox_bs, Ox, xrange(1):xrange(2));
        [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
        cx = cx(:); cy = cy(:);
        idx = I_mask(yrange(1):yrange(2),xrange(1):xrange(2));
        I_mask_t = zeros(sizeI);
        idxIn = cy >= 1 & cy <= sizeI(1) & cx >= 1 & cx <= sizeI(2);
        I_mask_t(sub2ind(sizeI, cy(idx(:)&idxIn), cx(idx(:)&idxIn))) = 1;
        set(hmask, 'AlphaData', I_mask_t*.2);
        %
        I = im2double(V(:,:,t));
        [E,~,Es] = fit_plane_energy(I_bs, I, Ox_bs, Ox, xrange, newyrange, dx, ...
                [xrange(1),leftpos(t,1); ddadX_init, ddadX; ddaeX_init, ddaeX; xrange(2), rightpos(t,1)]);
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



