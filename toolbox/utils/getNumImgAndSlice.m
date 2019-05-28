function [num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filepattern)
%GETNUMIMGANDSLICE get no. of frames and no. of slices
flist = dir(directory);
fidx = arrayfun(@(x)sscanf(x.name, filepattern), flist, 'UniformOutput', false);
fidx = fidx(cellfun(@length, fidx)==2);
maxidx = max([fidx{:}], [], 2);
if maxidx(1) > maxidx(2)
    switchTZ = false;
    num_img = maxidx(1);
    num_slice = maxidx(2);
else
    switchTZ = true;
    num_img = maxidx(2);
    num_slice = maxidx(1);
end
end

