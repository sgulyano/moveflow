function [ centers, radii, metric ] = detect_soma( I, opt )
%DETECT_SOMA use hough transform to estimate soma location

if nargin < 2;  opt = struct();  end;
if ~isfield(opt,'sigma');       opt.sigma     = 3;          end;
if ~isfield(opt,'adjval');      opt.adjval    = [0.1 .3];   end;
if ~isfield(opt,'windowgau');   opt.windowgau = [3 3];      end;
if ~isfield(opt,'radii');       opt.radii     = [15 30];    end;
if ~isfield(opt,'sensitiv');    opt.sensitiv  = 0.7;        end;
if ~isfield(opt,'edge_thr');    opt.edge_thr  = 0.3;        end;
if ~isfield(opt,'scale');       opt.scale     = 4;          end;


%%
Iadj = imadjust(I, opt.adjval, [0 1]);
h = fspecial('gaussian', opt.windowgau, opt.sigma);
f = imfilter(Iadj,h);
% f = f - min(f(:)); f = f ./ max(f(:));
% f = imdilate(imerode(f, ones(5)), ones(5));

fb = imresize(f, opt.scale);

[centers, radii, metric] = imfindcircles(fb, opt.radii, 'ObjectPolarity', 'bright', ...
        'Sensitivity', opt.sensitiv, ...
        'EdgeThreshold',opt.edge_thr, ...
        'Method','TwoStage');

centers = centers ./ opt.scale;
radii = radii ./ opt.scale;

% figure(2); imshow( f );
% viscircles(centers, radii,'EdgeColor','b');
% text(centers(:,1), centers(:,2),arrayfun(@num2str, metric, 'UniformOutput', false), 'Color', 'r', 'FontSize', 14) 
%%
% keyboard;
end

