function alldata = read_swc_file(filename)
%READ_SWC_FILE read neuron file(s) (swc file)
%   Input
%   filename - swc filename
%
%   Output
%   alldata - n-by-7 data in neuron file
%
if ~iscell(filename)
    filename = {filename};
end

alldata = [];
for i = 1:length(filename)
    disp(['Read SWC: ' filename{i}]);
    fid=fopen(filename{i}, 'r');
    if fid == -1 
        error('File could not be opened, check name or path.')
    end

    vnum = [];
    tline = '';
    while ischar(tline) && isempty(vnum)
        tline = fgetl(fid);
        vnum = sscanf(tline, '%d %d %f %f %f %f %d');
    end
    A = fscanf(fid,'%f');
    fclose(fid);
    data = [vnum'; reshape(A,7,length(A)/7)'];
    offset = size(alldata,1);
    data(data(:,1)>0,1) = data(data(:,1)>0,1) + offset;
    data(data(:,7)>0,7) = data(data(:,7)>0,7) + offset;
    alldata = [alldata; data];
end
end

