function soma = soma_manual_multi( V, t_init, soma_init, centers, radiis )
%SOMA_MANUAL manually change soma position
t = t_init;
n_neuron = size(soma_init,1);
soma = zeros([size(soma_init) size(V,3)]);
soma(:,:,t_init) = soma_init;
selected = 0;


f = figure('Position',[10 100 1080 720], 'NumberTitle','off', ...
        'MenuBar', 'none', 'Name','Tracking Soma');
    
uimenu(f,'Label','Help','Callback',['helpdlg({''Left click to select soma'',''w = up'',',...
        '''s = down'',''a = left'',''d = right'',''q/e = +/- radius'',''ctrl+r = remove'',',...
        '''space/right arrow = next frame'',''left arrow = previous frame'',',...
        '''Press Esc when finish''})']);
axes('Units','normalized','Position',[0.05,0.05,0.93,0.93],'FontSize',16);
himg = imshow(V(:,:,t), []); 
set(himg,'ButtonDownFcn',@mouseCallback);

title(['Tracking Soma at Frame ' num2str(t-1)]);
theta = -0.01:0.05:2*pi;
hold on
hsoma = cell(1,n_neuron);
for i = 1:n_neuron
    hsoma{i} = plot(soma_init(i,1)+soma_init(i,3)*cos(theta), soma_init(i,2)+soma_init(i,3)*sin(theta), 'g', 'LineWidth', 2);
end
hold off
set(f, 'WindowKeyPressFcn', @keyPressFcn);
uiwait(f);
pause(1);
try
    uiwait(msgbox('Done'));
catch
    warning('User hold spacebar. msgbox was closed before uiwait is activated.');
end
delete(f);


    function mouseCallback(hObject, eventdata)
        axesHandle  = get(hObject,'Parent');
        coordinates = get(axesHandle,'CurrentPoint'); 
        coordinates = coordinates(1,1:2);
        [ds, pos] = min( sqrt( sum((soma(:,1:2,t) - ones(n_neuron,1)*coordinates).^2,2) ) );
        if ds < max(soma(pos,3,t)*3, 20)
            if selected > 0
                set(hsoma{selected}, 'Color', 'g');
            end
            selected = pos;
            set(hsoma{selected}, 'Color', 'r');
        end
    end


    function keyPressFcn(hObject, eventdata)
        switch eventdata.Key
            case 'w'
                if selected > 0
                    soma(selected,2,t) = soma(selected,2,t) - 1;
                end
            case 's'
                if selected > 0
                    soma(selected,2,t) = soma(selected,2,t) + 1;
                end
            case 'a'
                if selected > 0
                    soma(selected,1,t) = soma(selected,1,t) - 1;
                end
            case 'd'
                if selected > 0
                    soma(selected,1,t) = soma(selected,1,t) + 1;
                end
            case 'q'
                if selected > 0
                    soma(selected,3,t) = soma(selected,3,t) - 1;
                end
            case 'e'
                if selected > 0
                    soma(selected,3,t) = soma(selected,3,t) + 1;
                end
            case 'r'
                if ismember('control', eventdata.Modifier) && selected > 0
                    soma(selected,:,t) = 0;
                end
            case {'space','rightarrow'}
                t = t + 1;
                if t > size(V,3) || t < 1, 
                    uiresume; 
                else
                    soma(:,:,t) = soma_match_multi( V(:,:,t), V(:,:,t-1), centers{t}, radiis{t}, soma(:,:,t-1));
%                     if flag1,
%                         uiresume(gcbf);
%                     end;
                    title(['Tracking Soma at Frame ' num2str(t-1)]);
                end
            case 'leftarrow'
                if t > 1
                    t = t - 1;
                    title(['Tracking Soma at Frame ' num2str(t-1)]);
                end
            case 'escape'
                uiresume(gcbf);
        end
        if t >= 1 && t <= size(V,3),
            himg.CData = V(:,:,t);
            for ii = 1:n_neuron
                hsoma{ii}.XData = soma(ii,1,t)+soma(ii,3,t)*cos(theta);
                hsoma{ii}.YData = soma(ii,2,t)+soma(ii,3,t)*sin(theta);
            end
            pause(0.05);
        end
    end
end

