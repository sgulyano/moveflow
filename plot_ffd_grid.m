function plot_ffd_grid( I1, I2, O1x, O1y, O2x, O2y )
%PLOT_FFD_GRID plot grid of FFD

subplot(1,2,1); imshow(I1,[]);
hold on
xx = [O1x nan(size(O1x,1),1)]'; yy = [O1y nan(size(O1x,1),1)]';
plot(xx(:), yy(:), 'g');
xx = [O1x; nan(1,size(O1x,2))]; yy = [O1y; nan(1,size(O1x,2))];
plot(xx(:), yy(:), 'g');
hold off

subplot(1,2,2); imshow(I2,[]);
hold on
xx = [O2x nan(size(O2x,1),1)]'; yy = [O2y nan(size(O2x,1),1)]';
plot(xx(:), yy(:), 'g');
xx = [O2x; nan(1,size(O2x,2))]; yy = [O2y; nan(1,size(O2x,2))];
plot(xx(:), yy(:), 'g');
hold off

end

