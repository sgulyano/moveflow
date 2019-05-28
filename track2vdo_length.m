function [ len, neu_len, len_baseline, time ] = track2vdo_length ...
        ( vdoname, Is, traces, t_bs, opt )
%TRACK2VDO save tracking results to video

if nargin < 4;  opt = struct();  end;
if ~isfield(opt,'pos');         opt.pos         = [0 500 800 500];  end;
if ~isfield(opt,'vec_sp');      opt.vec_sp      = 8;                end;
if ~isfield(opt,'soma_pos');    opt.soma_pos    = [];               end;
if ~isfield(opt,'soma_adj');    opt.soma_adj    = [];               end;
if ~isfield(opt,'tpf');         opt.tpf         = 58.2326 / 1000;   end;

col = [0 0 0 ;...
        0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        1.0000    0.0000    0.0000];
[num_img, num_trace] = size(traces);
num_img = num_img-1;
neu_len = zeros(size(traces));
time = (0:size(traces,1)-1) * opt.tpf;
len = zeros(size(traces));
len_baseline = zeros(size(traces));

outputVideo = VideoWriter(vdoname);
open(outputVideo);
h = figure(99); clf; set(gcf, 'Position', opt.pos);

len_bs = zeros(num_trace,1);
for j = 1:num_trace
    swc = traces{t_bs(j),j};
    prev = swc(:,7);
    idx = prev > 0; 
    n1 = swc(idx,1);
    n2 = prev(idx);
    len_bs(j) = sum(sqrt(sum((swc(n1,3:4) - swc(n2,3:4)).^2,2)));
end
% keyboard;
sizeI = [size(Is,1),size(Is,2)];
for t = 0:num_img
    I = Is(:,:,t+1);
    % get label volume
    label = zeros(sizeI);
    for num = 1:num_trace
        if ~isempty(traces{t+1,num})
            img = swc2pixel( traces{t+1,num}, sizeI );
            mask = imdilate(img, ones(3));
            
            swc = traces{t+1,num};
            prev = swc(:,7);
            idx = prev > 0; 
            n1 = swc(idx,1);
            n2 = prev(idx);
            neu_len(t+1,num) = sum(sqrt(sum((swc(n1,3:4) - swc(n2,3:4)).^2,2)));
            len_baseline(t+1,num) = len_bs(num);
            
            len(t+1,num) = abs(len_bs(num) - neu_len(t+1,num)) / neu_len(t+1,num);
            label(mask) = num;
        end
    end
    L = reshape(col(label+1,:), [sizeI 3]);
    mask = repmat(label==0, [1 1 3]);
    img = repmat(im2double(I), [1 1 3]);
    L(mask) = img(mask);
    subplot(2,1,1); imagesc((repmat(im2double(I), [1 1 3])+L)./2);
    if ~isempty(opt.soma_pos)
        for num = 1:num_trace
            viscircles(opt.soma_pos{num}(t+1,1:2), opt.soma_pos{num}(t+1,3), 'EdgeColor', col(num+1,:));
            
            if all(opt.soma_pos{num}(t+1,:) > 0)
                text(opt.soma_pos{num}(t+1,1)-5, opt.soma_pos{num}(t+1,2), num2str(num), ...
                        'Color','r','FontSize',16)
            end
        end
        for num = 1:num_trace
            viscircles(opt.soma_adj{num}(t+1,1:2), opt.soma_adj{num}(t+1,3), ...
                    'EnhanceVisibility', false,...
                    'EdgeColor', col(num+1,:), 'LineStyle', '--');
        end
    end
    title(sprintf('Time %.3f s, Frame %d', time(t+1), t)); drawnow;
    subplot(2,1,2);
    plot(time,len);
    xlabel('Time (s)');
    ylabel('|dD|/D');
    writeVideo(outputVideo, getframe(gcf));
end
close(outputVideo);
close(h);
disp(['Track Video saved to : ' vdoname]);
end

