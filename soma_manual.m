function [ cen, rad ] = soma_manual( img, cen, rad, init_t, inc, centers, radii )
%SOMA_MANUAL manually change soma position
t = init_t;
f = figure('Position',[100, 10, 1130, 580], 'NumberTitle','off', ...
        'MenuBar', 'none', 'Name','Tracking Soma');
    
uimenu(f,'Label','Help','Callback',['helpdlg({''w = up'',',...
        '''s = down'',''a = left'',''d = right'',''q/e = +/- radius'',',...
        '''space/right arrow = next frame'',''left arrow = previous frame'',',...
        '''Press Esc when finish''})']);
axes('Units','normalized','Position',[0.05,0.05,0.93,0.93],'FontSize',16);
himg = imshow(img(:,:,t), []); 

title(['Tracking Soma at Frame ' num2str(t-1)]);
theta = -0.01:0.05:2*pi;
hold on
hsoma = plot(cen(t,1)+rad(t)*cos(theta), cen(t,2)+rad(t)*sin(theta), 'r', 'LineWidth', 2);
hold off
set(f, 'WindowKeyPressFcn', @somaWindowKeyPressFcn);
uiwait(f);
pause(1);
try
    uiwait(msgbox('Done'));
catch
    warning('User hold spacebar. msgbox was closed before uiwait is activated.');
end
delete(f);


function somaWindowKeyPressFcn(hObject, eventdata, handles)
    switch eventdata.Key
        case 'w'
            cen(t,2) = cen(t,2) - 1;
        case 's'
            cen(t,2) = cen(t,2) + 1;
        case 'a'
            cen(t,1) = cen(t,1) - 1;
        case 'd'
            cen(t,1) = cen(t,1) + 1;
        case 'q'
            rad(t) = rad(t) - 1;
        case 'e'
            rad(t) = rad(t) + 1;
        case {'space','rightarrow'}
            t = t + inc;
            if t > size(img,3) || t < 1, 
                uiresume; 
            else
                [cen(t,:), rad(t), flag1] = soma_match(img(:,:,t), img(:,:,t-inc), ...
                        centers{t}, radii{t}, cen(t-inc,:), rad(t-inc));
                if flag1,
                    uiresume(gcbf);
                end;
                title(['Tracking Soma at Frame ' num2str(t-1)]);
            end
        case 'leftarrow'
            if t > 0
                t = t - inc;
                title(['Tracking Soma at Frame ' num2str(t-1)]);
            end
        case 'escape'
            uiresume(gcbf);
    end
    if t >= 1 && t <= size(img,3),
        himg.CData = img(:,:,t);
        hsoma.XData = cen(t,1)+rad(t)*cos(theta);
        hsoma.YData = cen(t,2)+rad(t)*sin(theta);
        pause(0.05);
    end
end
end

