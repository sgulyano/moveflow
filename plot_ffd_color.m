function plot_ffd_color( I, O1x, O1y, O2x, O2y, maxval )
%PLOT_FFD_COLOR Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    maxval = 5;
end

O_gradx = O2x - O1x;
O_grady = O2y - O1y;

dO_gradx = min(max(diff(O_gradx,[],2), -maxval), maxval) ./ maxval;
dO_grady = min(max(diff(O_grady,[],1), -maxval), maxval) ./ maxval;

%%
imshow(max(I,[],3),[]);
hold on
for i = 3:size(O1x,1)-1
    for j = 3:size(O1x,2)-1
        if O1y(i+1,j) <= size(I,1)
            if dO_grady(i,j) > 0
                plot(O1x(i:i+1,j), O1y(i:i+1,j,1), 'Color', [1-dO_grady(i,j) 1 0], 'LineWidth', 2);
            else
                plot(O1x(i:i+1,j), O1y(i:i+1,j,1), 'Color', [1 1+dO_grady(i,j) 0], 'LineWidth', 2);
            end
        end
        if O1x(i,j+1) <= size(I,2)
            if dO_gradx(i,j) > 0
                plot(O1x(i,j:j+1), O1y(i,j:j+1), 'Color', [1-dO_gradx(i,j) 1 0], 'LineWidth', 2);
            else
                plot(O1x(i,j:j+1), O1y(i,j:j+1), 'Color', [1 1+dO_gradx(i,j) 0], 'LineWidth', 2);
            end
        end
    end
end
hold off

end

