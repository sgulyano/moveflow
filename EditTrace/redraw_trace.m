function [ dendrite_plot ] = redraw_trace( swcdata, dendrite_plot )
%REDRAW_TRACE redraw trace plot on GUI

% remove old trace plot
for hdldendrite = dendrite_plot
    if hdldendrite == 0, continue, end
    delete(hdldendrite);
end
    
% redraw trace plot
hold on
nb = swcdata(:,[1 7]);
nb(any(nb < 1,2),:) = [];
G = sparse(nb(:,1), nb(:,2), 1, size(swcdata,1), size(swcdata,1));
G = G + G';
[cc, sizes] = components(G);
col = jet(max(5, length(sizes)+1));
dendrite_plot = zeros(1,length(sizes));
try
for j = 1:length(sizes)
    idx = all(cc(nb)==j,2);
    if ~any(idx), continue, end
    xx = [swcdata(nb(idx,1),3), swcdata(nb(idx,2),3), nan(sum(idx),1)]';
    yy = [swcdata(nb(idx,1),4), swcdata(nb(idx,2),4), nan(sum(idx),1)]';
    dendrite_plot(j) = plot(xx(:), yy(:), 'LineWidth', 1, 'Color', col(j,:));
end
%disp(dendrite_plot)
catch
    keyboard;
end
hold off

end

