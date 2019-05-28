function [ lines, direction ] = linegrowing( arr, tol )
%LINEGROWING Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2,
    tol = 3;
    nconsec = 3;
end
assert(isvector(arr));

%% detect line segments
nline = 2;
lines = [1 1];

sumval = -1;
st = 0;
numval = 0;
for i = 1:length(arr)
    if sumval < 0 || abs(sumval/numval - arr(i)) > tol
        if numval > nconsec
            lines(nline,:) = [st, st+numval-1];
            nline = nline + 1;
        end
        sumval = arr(i);
        numval = 1;
        st = i;
    else
        sumval = sumval + arr(i);
        numval = numval + 1;
    end
end

if numval > nconsec
	lines(nline,:) = [st, st+numval-1];
    nline = nline + 1;
end
lines(nline,:) = [length(arr) length(arr)];

%% merge adjacent line segments if they are too close
i = 1;
while i < size(lines,1)
    if lines(i,2)+nconsec >= lines(i+1,1)
        lines(i,2) = lines(i+1,2);
        lines(i+1,:) = [];
    else
        i = i + 1;
    end
end

%% compute mean
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

% %% sanity check
% figure, plot(arr, 'r.');
% hold on;
% plot(direction, 'm-');
% for i = 1:size(lines,1)
%     plot(lines(i,:), lines_mean(i) * ones(1,2), 'b-');
% end
% hold off;
% keyboard;

end

