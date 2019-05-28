function [ lines, direction ] = linegrowing( arr, opt )
%LINEGROWING find segments in arr that are stationary (values change within
%the tol). To avoid spurious segment, the number of consecutive stationary
%points (nconsec) must be met. The direction is determined by the
%difference between adjacent segments (-1 is left and 1 is right).
if nargin < 2;  opt = struct();  end;
if ~isfield(opt,'tol');             opt.tol     = 3;        end;
if ~isfield(opt,'nconsec');         opt.nconsec = 3;      	end;
if ~isfield(opt,'debug');           opt.debug   = false;    end;

% check input arr is a vector
assert(isvector(arr));

%% detect line segments
nline = 1;      % initialze with the empty array
lines = [1, 1];

sumval = -1;    % total sum of the segment (-1 for empty segment)
st = 0;         % starting index of the segment
numval = 0;     % number of values in the segment
for i = 1:length(arr)
    if sumval < 0 || abs(sumval/numval - arr(i)) > opt.tol
        if numval > opt.nconsec
            nline = nline + 1;
            lines(nline,:) = [st, st+numval-1];
        end
        sumval = arr(i);
        numval = 1;
        st = i;
    else
        sumval = sumval + arr(i);
        numval = numval + 1;
    end
end

if numval > opt.nconsec
    nline = nline + 1;
	lines(nline,:) = [st, st+numval-1];
end
lines(nline+1,:) = [length(arr), length(arr)];

%% merge adjacent line segments if they are too close
i = 1;
while i < size(lines,1) && size(lines,1) > 2
    if lines(i,2) + opt.nconsec >= lines(i+1,1)
        lines(i,2) = lines(i+1,2);
        lines(i+1,:) = [];
    else
        i = i + 1;
    end
end

%% compute mean of refined segments
lines_mean = zeros(size(lines,1),1);
for i = 1:size(lines,1)
    lines_mean(i) = mean(arr(lines(i,1):lines(i,2)));
end
lines_dir = sign(diff(lines_mean));

%% assign direction
direction = zeros(size(arr));
for i = 1:length(lines_dir)
    if i == 1
        direction(lines(i,1):lines(i+1,2)) = lines_dir(i);
    else
        direction(ceil(mean(lines(i,:))):lines(i+1,2)) = lines_dir(i);
    end
end

%% remove direction zero
[~,idx] = bwdist(direction);
direction = direction(idx);

%% sanity check
if opt.debug
    figure, plot(arr, 'r.');
    hold on;
    plot(direction, 'm-');
    for i = 1:size(lines,1)
        plot(lines(i,:), lines_mean(i) * ones(1,2), 'b-');
    end
    hold off;
    keyboard;
end
end

