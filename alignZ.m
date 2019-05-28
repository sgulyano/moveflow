function [zflips] = alignZ( V, swcs, track_time, soma_pos, configs, opt )
%ALIGNZAXIS GUI for adjusting Z axis of the neuron plane.

if nargin < 6;  opt = struct();  end;
if ~isfield(opt,'dx');          opt.dx = 4;                         end;
if ~isfield(opt,'pixelsize');   opt.pixelsize = 0.624;              end;
if ~isfield(opt,'tpf');         opt.tpf       = 58.2326 / 1000;     end;
if ~isfield(opt,'zflip');       opt.zflip     = [];                 end;

dx = opt.dx;
pixelsize = opt.pixelsize;

time = (0:size(soma_pos{1},1)-1) * opt.tpf;

Zmax = 30;%max(length(configs{track_time+1,1}), length(configs{track_time+1,2}))*dx*pixelsize;

midId = length(configs{track_time+1,1})+1;

if isempty(opt.zflip)
    zflips = cell(size(configs,1),1);
else
    zflips = opt.zflip;
end

%% init
ddad_somaX = round(swcs{1}(1,3));
ddae_somaX = round(swcs{2}(1,3));
somaX = round((ddad_somaX + ddae_somaX)/2);
[Ox_bs, Oz_bs] = ffd_init_from_config(configs(track_time+1,:), somaX, dx);
Ox_bs = Ox_bs .* pixelsize; Oz_bs = Oz_bs .* pixelsize;
if isempty(opt.zflip)
    zflips{track_time+1} = ones(size(Oz_bs));
end

% init plot
f = figure(1);
subplot(2,1,1); hxz = plot(Ox_bs, Oz_bs, '-or', 'MarkerSize',12, 'MarkerFaceColor', 'r');
hold on;
hmid = scatter(Ox_bs(midId), Oz_bs(midId), 200, 'bo', 'filled');
hcon = plot(Ox_bs([1 end]), -9.56*ones(1,2), '--b');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('X (\mum)');
ylabel('Z (\mum)');
ylim([-Zmax Zmax]);
htitle = title(sprintf('Time %.3f s, Frame %d', time(track_time+1), track_time));

subplot(2,1,2); hxz_old = plot(Ox_bs, Oz_bs, '-or', 'MarkerSize',12, 'MarkerFaceColor', 'r');
hold on;
hmid_old = scatter(Ox_bs(midId), Oz_bs(midId), 200, 'bo', 'filled');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('X (\mum)');
ylabel('Z (\mum)');
ylim([-Zmax Zmax]);

% setup GUI
set(f,'ButtonDownFcn', @(~,~)disp('figure'), 'HitTest', 'off')
set(f,'WindowKeyPressFcn', @keyboardCallback);
set(hxz,'ButtonDownFcn', @mouseCallback, 'PickableParts', 'all', 'HitTest', 'on');

% adjust Z-axis in each frame
for t = track_time+2:size(V,3)
    if isempty(configs{t,1})
        break;
    end
    if isempty(opt.zflip)
        zflips(t) = zflips(t-1);
    end
    ddadX = soma_pos{1}(t,1);
    ddaeX = soma_pos{2}(t,1);
    newsomaX = (ddadX + ddaeX)/2;
    % plot the original plane
    [Ox_old, Oz_old] = ffd_init_from_config(configs(t,:), newsomaX, dx);
    Ox_old = Ox_old .* pixelsize; Oz_old = Oz_old .* pixelsize;
    set(hxz_old, 'XData', Ox_old);
    set(hxz_old, 'YData', Oz_old);
    set(hmid_old, 'XData', Ox_old(midId));
    set(hmid_old, 'YData', Oz_old(midId));
    set(htitle, 'String', sprintf('Time %.3f s, Frame %d', time(t), t-1));
    
    % plot the new plane
    [Ox, Oz] = ffd_init_from_config(configs(t,:), newsomaX, dx, zflips{t});
    Ox = Ox .* pixelsize; Oz = Oz .* pixelsize;
    set(hcon, 'XData', Ox([1 end]));
    set(hxz, 'XData', Ox);
    set(hxz, 'YData', Oz);
    set(hmid, 'XData', Ox(midId));
    set(hmid, 'YData', Oz(midId));
    
    if all(Oz == 0)
        pause(0.1);
    else
        uiwait(gcf);
    end
end

	function mouseCallback(src, evt)
        figHandle = ancestor(src, 'figure');
        clickType = get(figHandle, 'SelectionType');

        if strcmp(clickType, 'alt')
            disp('right click action goes here!');
            uiresume(gcbf)
        else
            axesHandle  = get(src,'Parent');
            coordinates = get(axesHandle,'CurrentPoint'); 
            coordinates = coordinates(1,1:2);
            disp(coordinates);
            
            % find selected point
            ds = (Ox - coordinates(1)).^2 + (Oz - coordinates(2)).^2;
            [~, pos] = min(ds);
            zflips{t}(pos) = zflips{t}(pos) * -1;
            
            % update after flip z
            [Ox, Oz] = ffd_init_from_config(configs(t,:), newsomaX, dx, zflips{t});
            Ox = Ox .* pixelsize; Oz = Oz .* pixelsize;
            set(hxz, 'XData', Ox);
            set(hxz, 'YData', Oz);
            set(hmid, 'XData', Ox(midId));
            set(hmid, 'YData', Oz(midId));
        end
    end

    function keyboardCallback(hObject, eventdata, handles)
        switch eventdata.Key
            case 'space'
                uiresume(gcbf)
        end
    end

end

