function imgflow2vdo( vdoname, I, U, V, opt )
%IMGFLOW2VDO save image stack along with vector flow in video
% Input:
%   I - image stack
%   U - X-direction of vector field
%   V - Y-direction of vector field

if nargin < 4
    U = [];
    V = [];
end

if nargin < 5;  opt = struct();  end;
if ~isfield(opt,'pos');     opt.pos     = [0 500 1000 400];    end;
if ~isfield(opt,'vec_sp');  opt.vec_sp  = 8;    end;

FNUM = size(I,3);
U = padarray(U, [0 0 max(FNUM-size(U,3),0)], 0, 'post');
V = padarray(V, [0 0 max(FNUM-size(V,3),0)], 0, 'post');

outputVideo = VideoWriter(vdoname, 'MPEG-4');
open(outputVideo);
h = figure(99); clf; set(gcf, 'Position', opt.pos);

vecmask = false(size(I,1), size(I,2));
vecmask(opt.vec_sp:opt.vec_sp:size(I,1),opt.vec_sp:opt.vec_sp:size(I,2)) = true;
for t = 1:FNUM
    if isempty(U)
        if ~exist('h1', 'var')
            h1 = imagesc(I(:,:,t)); colormap gray;
        else
            set(h1, 'CData', I(:,:,t));
        end
        title(t-1);
        axis off;
        drawnow;
    else
        im2 = repmat(I(:,:,t), [1 1 3]);
        sc = max(max(sqrt(U(:,:,t).^2+V(:,:,t).^2)));

        if ~exist('h1', 'var') || ~ishandle(h1)
            [h1, h2] = plotFlow(U(:,:,t), V(:,:,t), im2, opt.vec_sp, sc/2);
        else
            set(h1, 'CData', im2);
            Utmp = U(:,:,t); Utmp(~vecmask) = 0;
            set(h2, 'UData', Utmp);
            Vtmp = V(:,:,t); Vtmp(~vecmask) = 0;
            set(h2, 'VData', Vtmp);
            set(h2, 'AutoScaleFactor', sc/2);
        end
        title(t-1);
        axis off
    end
    drawnow;
    writeVideo(outputVideo, getframe(gcf));
end
close(outputVideo);
close(h);
disp(['Video saved to : ' vdoname]);

end

