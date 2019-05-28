function main()
%% User Param
DEBUG = false;
DEPLOY = true;
addpath(genpath('toolbox/utils'))
addpath(genpath('fix_slice'))
if ~DEPLOY
    addpath(genpath('toolbox/utils'))
    addpath(genpath('fix_slice'))
    % run before deploy
    txtfolder_default = '/Users/sarun/Desktop/LipingsData/20180404mcdGFP_Curvature/20180404-L6-B';
    txtfilename_default = '20180404-L6-B_z%02d_t%02d.tif';
    txtsave_default = fullfile(pwd, 'test');
    txtswcddad_default = '/Users/sarun/Desktop/moveflow/test';
    txtswcddae_default = '/Use             rs/sarun/Desktop/moveflow/test';
    save('config.mat','txtfolder_default','txtfilename_default','txtsave_default',...
            'txtswcddad_default','txtswcddae_default');
end

m = matfile('config.mat','Writable',true);

%% =====================================
%  Open dialog for select image stacks
%  =====================================
f = figure('Visible','on','Position',[360,500,600,200], ...
        'MenuBar','none','ToolBar','none',...
        'Name','Select folder and filename','NumberTitle','off');

% Construct the components.
% input directory selector
uicontrol('Style','text',...
        'Position',[15 160 80 25],...
        'TooltipString', ['The directory containing all slices in the image stack.' char(10) ...
        'Each image file is one slice; and one directory per stack.'],...
        'String','Directory','FontWeight','Bold','FontSize',14);
txtfolder = uicontrol('Style','text',...
        'Position',[100 130 400 60],'HorizontalAlignment','Left',...
        'String',m.txtfolder_default,'FontSize',12);
uicontrol('Style','pushbutton',...
        'String','Open','Position',[525,160,60,25],'FontSize',14,...
        'Callback',@openCallback);

% editbox for pattern of input filename 
uicontrol('Style','text',...
        'Position',[15 100 80 25],...
        'TooltipString', ['The pattern of image filenames. It is the filename, where slice number' char(10) ...
        'is replaced by %d. If the slice number has leading zeroes, then ' char(10) ...
        'put 0 followed by the number of digit like z%03d.tif for z001.tif'],...
        'String','Filename','FontWeight','Bold','FontSize',14);
txtfilename = uicontrol('Style','edit',...
        'Position',[100 100 350 25],'HorizontalAlignment','Left',...
        'String',m.txtfilename_default,'FontSize',12);

% save location selector
uicontrol('Style','text',...
        'Position',[15 50 80 40],...
        'TooltipString', ['The program will produce 3 outputs with the same name ' ...
        'but different extensions (.mat, .swc, and .log). ' char(10) ...
        'The most important one is SWC file, which you will need for editing and further analysis.'],...
        'String','Save Directory','FontWeight','Bold','FontSize',14);
txtsave = uicontrol('Style','text',...
        'Position',[100 35 400 60],'HorizontalAlignment','Left',...
        'String',m.txtsave_default,'FontSize',12);
uicontrol('Style','pushbutton',...
        'String','Open','Position',[525,65,60,25],'FontSize',14,...
        'Callback',@saveCallback);

% align stack button
uicontrol('Style','pushbutton',...
        'String','Align Stack','Position',[50,10,160,30],'FontSize',16,...
        'Callback',@alignstackCallback);

% pick frame button
uicontrol('Style','pushbutton',...
        'String','Pick Frame','Position',[240,10,160,30],'FontSize',16,...
        'Callback',@pickframeCallback);

% next button
uicontrol('Style','pushbutton',...
        'String','Next','Position',[430,5,100,40],'FontSize',16,...
        'Callback',@nextCallback);

    %% =====================================
    %  Callback functions
    %  =====================================
    function openCallback(~, ~)
        directory = uigetdir(m.txtfolder_default, 'Select Directory of Image Sequence');
        if ischar(directory)
            set(txtfolder, 'String', directory);
            m.txtfolder_default = directory;
        end
    end

    function saveCallback(~, ~)
        selpath = uigetdir(m.txtsave_default, 'Select Save Directory');
        if ischar(selpath)
            set(txtsave, 'String', selpath);
            m.txtsave_default = selpath;
        end
    end

    function alignstackCallback(~,~)
        directory = txtfolder.String;
        filename = txtfilename.String;
        fname = fullfile(directory, sprintf(filename,1,1));
        if exist(fname, 'file') ~= 2
            errordlg(['Cannot find the selected file. (' fname ')']);
            return
        end
        m.txtfilename_default = txtfilename.String;
        
        options.Default = 'Cancel';
        options.Interpreter = 'tex';
        answer = questdlg(['\fontsize{14} Make sure you pick the directory and type filename pattern of '...
                'the original image stack sequence. The aligned sequence will be stored in '...
                'the Save Directory and existing files may be replaced. Do you want to proceed?'], ...
                'Align Image Stack Sequence', ...
                'Proceed','Cancel',options);
        if strcmp(answer, 'Proceed')
            % copy specific frame for manually tracing the initial neurite model
            savedir = txtsave.String;
            aligndir = fullfile(savedir, 'aligned');
            if ~exist(aligndir, 'dir')
                mkdir(aligndir);
            end
            [num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);
            fix_slice_func(directory, filename, aligndir, num_img, num_slice, switchTZ);
            set(txtfolder, 'String', aligndir);
        end
    end

    function pickframeCallback(~, ~)
        % copy specific frame for manually tracing the initial neurite model
        directory = txtfolder.String;
        filename = txtfilename.String;
        savedir = txtsave.String;

        m.txtfilename_default = txtfilename.String;

        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end

        % let user pick the reference frame number
        prompt = {'Enter the reference frame number'};
        title = 'Input';
        dims = [1 35];
        definput = {'1'};
        answer = inputdlg(prompt,title,dims,definput);

        % validate input reference frame number: is number or not?
        if isempty(answer)
            errordlg('No reference frame is given');
            return;
        else
            [num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);
            t = str2double(answer{1});
        end

        % validate input reference frame number: is in a valid range of value?
        if isnan(t) || t < 0 || t > num_img
            errordlg('Invalid reference frame is given');
            return;
        else
            % create a new folder to store the frame image stack
            pickframefile = fullfile(savedir, ['ref_t' num2str(t) '.tif']);
            if exist(pickframefile, 'file') == 2
                answer = questdlg([pickframefile ' already exists. Do you want to replace it?'], ...
                        'File exists', ...
                        'Yes','No','No');
                if strcmp(answer, 'No')
                    return;
                end
            end

            % copy the specific frame
            disp('Write reference frame to file.');
            for z = 0:num_slice
                if switchTZ
                    sourcefile = fullfile(directory, sprintf(filename,z,t));
                else
                    sourcefile = fullfile(directory, sprintf(filename,t,z));
                end
                X = imread(sourcefile);
                if z==0
                    imwrite(X,pickframefile)
                else
                    imwrite(X,pickframefile,'WriteMode','append')
                end
            end
            % done
            msgbox(['The reference frame ' num2str(t) ' is written in ' pickframefile '.']);
        end
    end


    function nextCallback(~, ~)
    % ask users about the image sequence before tracking can be done
    directory = txtfolder.String;
    filename = txtfilename.String;
    savedir = txtsave.String;
    savename = strsplit(filename, '_');
    savename = savename{1};
    
    m.txtfilename_default = txtfilename.String;
    
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
    
    %% create the dialog asking for input
    d = dialog('Position',[300 300 350 250],'Name','Input');
    % let user pick the reference frame number
    uicontrol('Parent',d,...
            'Style','text',...
            'Position',[10 200 250 40],...
            'FontSize', 14,...
            'HorizontalAlignment','Left',...
            'FontWeight','Bold',...
            'String','Number of ddaD and ddaE pairs');
    txtnumpair = uicontrol('Parent',d,...
            'Style','edit',...
            'Position',[10 190 70 25],...
            'FontSize', 14,...
            'String','1');
    % ask the movement direction of the image seqeuence
    uicontrol('Parent',d,...
            'Style','text',...
            'Position',[10 140 250 40],...
            'FontSize', 14,...
            'HorizontalAlignment','Left',...
            'FontWeight','Bold',...
            'String','Movement direction');
    bg = uibuttongroup('Parent',d,...
            'BorderWidth', 0,...
            'Visible','on',...
            'Units','pixels',...
            'Position',[10 120 350 40],...
            'SelectionChangedFcn','');
    rb = uicontrol(bg,'Style','radiobutton',...
            'String','Backward',...
            'FontSize', 14,...
            'Position',[10 10 140 30],...
            'HandleVisibility','off');
    rf = uicontrol(bg,'Style','radiobutton',...
            'String','Forward',...
            'FontSize', 14,...
            'Position',[150 10 140 30],...
            'HandleVisibility','off');
    % next button
    uicontrol('Parent',d,...
           'Position',[140 10 70 25],...
           'String','Next',...
           'FontSize', 14,...
           'Callback',@nextbtnCallback);
       
        %% callback function
        function nextbtnCallback(~,~)
            num_pairs = str2double(txtnumpair.String);
            if isnan(num_pairs) || num_pairs <= 0
                errordlg('Invalid number of ddaD/ddaE pairs');
                return;
            else
                if rb.Value
                    forwardmovement = false;
                else
                    forwardmovement = true;
                end
            end
            delete(d);

            %% =====================================
            %  Main Flow of GUI Start HERE
            %  =====================================
            % ask for SWC file of initial neurite model for each frame
            FILE_LOAD = false;
            if exist(fullfile(savedir, 'user_input.mat'), 'file') == 2
                answer = questdlg([fullfile(savedir, 'user_input.mat') ' already exists. Do you want to use the previous setting?'], ...
                        'Previous setting exists', ...
                        'Yes','No','Yes');
                if strcmp(answer, 'Yes')
                    load(fullfile(savedir, 'user_input.mat'), ...
                            'params', 'soma_cen', 'soma_rad', 'forwardmovement', ...
                            'initswctime', 'ddad_swc', 'ddae_swc', 'dd_plane');
                    FILE_LOAD = true;
                end
            end
            
            % read the middle frame for tuning parameter
            [num_img, num_slice, switchTZ] = getNumImgAndSlice(directory, filename);
            t = floor(num_img / 2);
            
            if ~FILE_LOAD
                % get neuron pairs
                initswctime = zeros(1,num_pairs);
                ddad_swc = cell(1,num_pairs);
                ddae_swc = cell(1,num_pairs);
                for i = 1:num_pairs
                    [initswctime(i), ddad_swc{i}, ddae_swc{i}] = getNeuronPair(i, m);
                end
            
                % read image stack for tuning parameters
                gfpinfo = imfinfo(fullfile(directory, sprintf(filename,0,0)));
                disp('Tuning parameter');
                V = zeros(gfpinfo.Height, gfpinfo.Width, num_slice+1, 'uint8');
                for z = 0:num_slice
                    if switchTZ
                        X = imread( fullfile(directory, sprintf(filename,z,t)) );
                    else
                        X = imread( fullfile(directory, sprintf(filename,t,z)) );
                    end
                    X = single(X);
                    X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));

                    n = 3;
                    Idouble = im2double(X);
                    avg = mean2(Idouble);
                    sigma = std2(Idouble);
                    X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);

                    V(:,:,z+1) = X;
                end

                % tuning parameters
                params = tuneParam( double(V) );
            end
            [gfp, gfp1sl] = readImg(directory, filename, params);
            if isempty(gfp)
                errordlg(['Cannot find the selected file. (' sprintf(filepattern, 1, 1) ')']);
                return
            end
            delete(f);
            % track neuron
            if ~FILE_LOAD
                [soma_cen, soma_rad, dd_plane] = neuron_tracking(savedir, gfp, gfp1sl, num_pairs, initswctime, ddad_swc, ddae_swc);
                
                if DEBUG
                    load(fullfile(savedir, 'user_input.mat'))
                else
                    save(fullfile(savedir, 'user_input.mat'), 'directory', 'filename', ...
                            'params', 'soma_cen', 'soma_rad', 'forwardmovement', ...
                            'initswctime', 'ddad_swc', 'ddae_swc', 'dd_plane');
                end
            else
                % view tracking results
                for i = 1:num_pairs
                    plot_track( gfp, dd_plane(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)] );
                end
            end
            
            ddad_lefts = zeros(1,num_pairs);
            for i = 1:num_pairs
                ddad_lefts(i) = sign(ddad_swc{i}(1,3) - ddae_swc{i}(1,3));
            end
            ddaD_Left = mode(ddad_lefts) == -1;
            
            % convert to excel            
            A = load('thetain.mat', 'thetain');
            opt.is_gui = true;
            opt.tval = A.thetain;
            opt.ddad_left = ddaD_Left;
            if forwardmovement
                A = load('ddaeforward_back.mat', 'ddaeforward_back_sem', 'ddaeforward_back_mean');
                opt.tsem = A.ddaeforward_back_sem;
                opt.tmean = A.ddaeforward_back_mean;
            else
                A = load('ddadbackward_back.mat', 'ddadbackward_back_sem', 'ddadbackward_back_mean');
                opt.tsem = A.ddadbackward_back_sem;
                opt.tmean = A.ddadbackward_back_mean;
            end
            [~, directions] = cellfun(@(x)linegrowing(x(:,1)), soma_cen, 'UniformOutput', false);
            opt.directions = directions;
            convert2excel(savedir, savename, gfp, num_img, initswctime, ddad_swc, ddae_swc, soma_cen, soma_rad, dd_plane, opt);
            if DEBUG
                keyboard;
            end
        end
    end
end

%% =====================================
%  Tracking function
%  ===============        ======================
function [soma_cen, soma_rad, results] = neuron_tracking(savedir, gfp, gfp1sl, num_pairs, initswctime, ddad_swc, ddae_swc)
%% track soma
[soma_cen, soma_rad] = track_soma(gfp1sl, initswctime, [ddad_swc; ddae_swc]);

%% track a pair of ddaD and ddaE by fitting a plane through folding
results = struct([]);
for i = 1:num_pairs
     [configs, leftpos, rightpos, spline] = ffd_gui_both(gfp, {ddad_swc{i}, ddae_swc{i}}, ...
            initswctime(i), {soma_cen{1,i}, soma_cen{2,i}});
%     videoname = fullfile(savedir, sprintf('model_both%d.mp4', i));
%     ffd_plot_both(videoname, gfp, {ddad_swc{i}, ddae_swc{i}}, ...
%             initswctime(i), {soma_cen{1,i}, soma_cen{2,i}}, ...
%             configs, spline, struct('is_gui', true));
    results(i).configs = configs;
    results(i).leftpos = leftpos;
    results(i).rightpos = rightpos;
    results(i).spline = spline;
    
    plot_track( gfp, results(i), ddad_swc{i}, ddae_swc{i}, soma_cen(:,i), ['Track Result Pair No. ' num2str(i)] );
    
    answer = questdlg('Do you want to save this result?', ...
                'Save result', ...
                'Yes','No' ,'Yes  ');
    if strcmp(answer, 'No')
        continue;
    end
end
end


%% ============================================
%  Dialog for getting SWC file of neuron pair
%  ============================================
function [frame_num, ddad_swc, ddae_swc] = getNeuronPair(i, m)
%%
d = dialog('Position',[300 300 350 300],'Name',['Neuron Pair ' num2str(i)]);
uicontrol('Parent',d,...
        'Style','text',...
        'Position',[10 250 250 40],...
        'FontSize', 14,...
        'HorizontalAlignment','Left',...
        'FontWeight','Bold',...
        'String','Frame Number');
txtframenum = uicontrol('Parent',d,...
        'Style','edit',...
        'Position',[10 240 70 25],...
        'FontSize', 14,...
        'String','1');

uicontrol('Parent',d,...
        'Style','text',...
        'Position',[10 205 250 25],...
        'FontSize', 14,...
        'HorizontalAlignment','Left',...
        'FontWeight','Bold',...
        'String','SWC File of ddaD');
uicontrol('Parent',d,...
        'Position',[230 205 70 25],...
        'String','Open',...
        'FontSize', 14,...
        'Callback',@openddadCallback);
txtswcddad = uicontrol('Parent',d,...
        'Style','text',...
        'HorizontalAlignment','Left',...
        'String', fullfile(m.txtswcddad_default, 't01_ddaD1.swc'),...
        'Position',[10 140 330 65],...
        'FontSize', 12);

uicontrol('Parent',d,...
        'Style','text',...
        'Position',[10 115 250 25],...      
        'FontSize', 14,...
        'HorizontalAlignment','Left',...
        'FontWeight','Bold',...
        'String','SWC File of ddaE');
uicontrol('Parent',d,...
        'Position',[230 115 70 25],...
        'String','Open',...
        'FontSize', 14,...          
        'Callback',@openddaECallback);
txtswcddae = uicontrol('Parent',d,...
        'Style','text',...
        'HorizontalAlignment','Left',...
        'String',fullfile(m.txtswcddae_default, 't01_ddaE1.swc'),...
        'Position',[10 50 330 65],...
        'FontSize', 12);

uicontrol('Parent',d,...
        'Position',[140 10 70 25],...
        'String','Okay',...
        'FontSize', 14,...
        'Callback',@okayCallback);
while true
    uiwait(d);
    frame_num = str2double(txtframenum.String);
    if isnan(frame_num)
        errordlg('Frame Number is invalid');
    end
    if exist(txtswcddad.String, 'file') ~= 2
        errordlg([txtswcddad.String ' does not exist.']);
        continue
    end
    ddad_swc = read_swc_file(txtswcddad.String);
    if exist(txtswcddae.String, 'file') ~= 2
        errordlg([txtswcddae.String ' does not exist.']);
        continue
    end
    ddae_swc = read_swc_file(txtswcddae.String);
    close(d);
    break
end

    function openddadCallback(~,~)
        [file, path] = uigetfile(fullfile(m.txtswcddad_default, '*.swc'),...
                'SWC File Selector for ddaD');
        if ischar(file) && ischar(path)
            if exist([path file], 'file') == 2
                set(txtswcddad, 'String', [path file]);
                m.txtswcddad_default = path;
            else
                errordlg(['Cannot find the selected file. (' path file ')']);
            end
        end
    end

    function openddaECallback(~,~)
        [file, path] = uigetfile(fullfile(m.txtswcddae_default, '*.swc'),...
                'SWC File Selector for ddaE');
        if ischar(file) && ischar(path)
            if exist([path file], 'file') == 2
                set(txtswcddae, 'String', [path file]);
                m.txtswcddae_default = path;
            else
                errordlg(['Cannot find the selected file. (' path file ')']);
            end
        end
    end
    
    function okayCallback(~,~)
        % validate input
        frame_num = str2double(txtframenum.String);
        if isnan(frame_num) || frame_num < 0
            errordlg('Invalid frame number');
        else
            uiresume(gcbf);
        end
    end
end


%% =====================================
%  GUI for tracking soma manually
%  =====================================
function [soma_cen, soma_rad] = track_soma(gfp1sl, initswctime, swcs)
type_str = {'ddaD', 'ddaE'};
%% detect soma candidates using hough transform
num_img = size(gfp1sl,3);
disp('Tracking soma... started');
soma_opt = struct('radii',[5 8], 'sensitiv',.95, 'edge_thr',.1);
centers = cell(num_img+1,1);
radii = cell(num_img+1,1);
    
fwait = waitbar(0,'Finding soma candidates in every frame');
for t = 1:num_img
    [centers{t}, radii{t}] = imfindcircles(gfp1sl(:,:,t), soma_opt.radii, ...
            'ObjectPolarity', 'bright', ...
            'Sensitivity', soma_opt.sensitiv, ...
            'EdgeThreshold',soma_opt.edge_thr, ...
            'Method','TwoStage');
    waitbar(t/num_img,fwait,['Finding soma candidates in frame ' num2str(t) '/' num2str(num_img)]);
end
close(fwait);

%% track soma for each trace
soma_cen = cell(size(swcs));
soma_rad = cell(size(swcs));
for i = 1:numel(swcs)
    num = floor((i-1)/size(swcs,1))+1;
    track_time = initswctime(num);
    pnt = swcs{i}(:,3:5)+1;

    soma_cen{i} = zeros(num_img, 2); 
    soma_rad{i} = zeros(num_img, 1);
    
    [~, pos] = min(sum(bsxfun(@minus, centers{track_time+1}, pnt(1,1:2)).^2,2));
    soma_cen{i}(track_time+1,:) = centers{track_time+1}(pos,:);
    soma_rad{i}(track_time+1) = radii{track_time+1}(pos);
    
    f = figure('Position', [100, 10, 1130, 580], 'NumberTitle','off', ...
            'Name',['Initial soma position of ' type_str{mod(i-1,2)+1} ' no. ' ...
            num2str(num) ' at frame ' num2str(track_time)]);
    uicontrol('Style','pushbutton','String','Start','FontSize',16,...
            'Units','normalized','Position',[0.01,0.525,0.05,0.1],...
            'Callback','uiresume(gcbf)');
    axes('Units','normalized','Position',[0.1,0.05,0.88,0.93],'FontSize',16);

    EV_plot_img(gfp1sl(:,:,track_time+1), swcs{i}); 
    hold on;
    EV_plot_img([], swcs{i});

    title('Red circle is the initial soma, blue ones are candidates, and cyan lines is the initial trace.');
    viscircles(centers{track_time+1}, radii{track_time+1}, 'EdgeColor', 'b');
    viscircles(soma_cen{i}(track_time+1,:), soma_rad{i}(track_time+1), 'EdgeColor', 'r');
    drawnow;
    hold off;
    uiwait(f);
    close(f);
    % track forward
    [soma_cen{i}, soma_rad{i}] = soma_manual( gfp1sl, soma_cen{i}, soma_rad{i}, track_time+1, 1, centers, radii );
end
end