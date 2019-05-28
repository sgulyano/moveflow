function [ activity_sum, activity_mean, activity_count ] = track2vdo( vdoname, Is, traces, opt )
%TRACK2VDO save tracking results to video

if nargin < 4;  opt = struct();  end;
if ~isfield(opt,'pos');         opt.pos         = [0 500 800 500];  end;
if ~isfield(opt,'vec_sp');      opt.vec_sp      = 8;                end;
if ~isfield(opt,'soma_pos');    opt.soma_pos    = [];               end;

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
activity_sum = zeros(size(traces));
activity_mean = zeros(size(traces));
activity_count = zeros(size(traces));

outputVideo = VideoWriter(vdoname);
open(outputVideo);
h = figure(99); clf; set(gcf, 'Position', opt.pos);

sizeI = [size(Is,1),size(Is,2)];
for t = 0:num_img
    I = Is(:,:,t+1);
    % get label volume
    label = zeros(sizeI);
    for num = 1:num_trace
        if ~isempty(traces{t+1,num})
            img = swc2pixel( traces{t+1,num}, sizeI );
            mask = imdilate(img, ones(3));
            activity_sum(t+1,num) = sum(double(I(mask)));
            activity_count(t+1,num) = sum(mask(:));
            activity_mean(t+1,num) = mean(double(I(mask)));
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
        end
    end
    title(num2str(t)); drawnow;
    subplot(2,1,2);
    plot(activity_sum);
    xlabel('Time Step');
    ylabel('Activity');
    writeVideo(outputVideo, getframe(gcf));
end
close(outputVideo);
close(h);
disp(['Track Video saved to : ' vdoname]);

end

