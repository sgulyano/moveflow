function [ cen, rad ] = soma_manual( img, cen, rad, init_t, inc, centers, radii )
%SOMA_MANUAL manually change soma position
t = init_t;
f = figure(5); himg = imshow(img(:,:,t), []); title(t-1);
set(gcf, 'pos', [10 100 1400 600]);
theta = -0.01:0.05:2*pi;
hold on
hsoma = plot(cen(t,1)+rad(t)*cos(theta), cen(t,2)+rad(t)*sin(theta), 'r', 'LineWidth', 2);
hold off
set(f, 'WindowKeyPressFcn', @somaWindowKeyPressFcn);
uiwait;
delete(f);
pause;
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
                    uiresume;
                end;
                title(t-1);
            end
        case 'leftarrow'
            t = t - inc;
            title(t-1);        
        case 'escape'
            uiresume;
    end
    if t >= 1 && t <= size(img,3),
        himg.CData = img(:,:,t);
        hsoma.XData = cen(t,1)+rad(t)*cos(theta);
        hsoma.YData = cen(t,2)+rad(t)*sin(theta);
        pause(0.05);
    end
%     disp(eventdata.Key);
end
end

