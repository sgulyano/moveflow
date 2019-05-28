addpath(genpath('../toolbox/utils'))
addpath(genpath('nonrigid_version23'));
addpath('..')

base_dir = 'C:\Users\sguly\Downloads\2018 New data for dendrites\2018 New data for dendrites\20181115UASmCD8GFP;221Gal4,tmc1';
save_dir = 'C:\Users\sguly\Documents\20181115UASmCD8GFP;221Gal4,tmc1_align';

flist = dir(base_dir);
for f = flist'
    % get zip file
    [~, fn, fext] = fileparts(f.name);
    if strcmp(fext, '.zip')
        disp(['---- ' f.name ' ----'])
        input_zip = fullfile(base_dir, f.name);
        output_dir = fullfile(base_dir, fn);
        % unzip the sequence zip file
        fns = unzip(input_zip, output_dir);
        
        [directory, fnsn, fnsext] = fileparts(fns{1});
        
        z_str = regexp(fnsn, '_z\d+', 'match');
        z_str = z_str{1};
        z_val = length(z_str(3:end));
        
        t_str = regexp(fnsn, '_t\d+', 'match');
        t_str = t_str{1};
        t_val = length(t_str(3:end));
        % get file pattern and number of images/slices
        pattern = strrep(fnsn, z_str, ['_z%0' num2str(z_val) 'd']);
        pattern = strrep(pattern, t_str, ['_t%0' num2str(t_val) 'd']);
        filename = [pattern, fnsext];
        [num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);
        
        % align image stacks
        aligndir = fullfile(save_dir, [fn '_aligned']);
        if ~exist(aligndir, 'dir')
            mkdir(aligndir);
            fix_slice_func(directory, filename, aligndir, num_img, num_slice, switchTZ);
            disp('Done')
        end
        
        % zip files back for uploading
        if ~exist([aligndir '.zip'], 'file')
            zip([aligndir '.zip'], aligndir);
            disp('Zipped Done')
        end
    end
end

