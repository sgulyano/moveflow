function [ params ] = tuneParam( V )
%TUNEPARAM let user tunes parameters
Vproc = im2double(V);
params = [];

%% define interface
f = figure('Visible','on','Position',[100, 10, 1080, 580], ...
        'Name','Tune Parameter','NumberTitle','off');

axes('Units','normalized','Position',[0.3,0.15,0.65,0.75],'FontSize',16);
himg = imshow(max(Vproc,[],3),[]);

hview = uibuttongroup('Visible','on',...
        'Units','normalized','Position',[0.48 0.92 0.3 0.05],...
        'SelectionChangedFcn',@viewSelectCallback);
              
% Create three radio buttons in the button group.
h2d = uicontrol(hview,'Style','radiobutton',...
        'String','2D','FontSize',14,...
        'Units','normalized','Position',[0.05 0.05 0.4 0.9],...
        'TooltipString', 'Show intensity maximum projection',...
        'HandleVisibility','off');
              
h3d = uicontrol(hview,'Style','radiobutton',...
        'String','3D','FontSize',14,...
        'Units','normalized','Position',[0.55 0.05 0.4 0.9],...
        'TooltipString', 'Show each slice',...
        'HandleVisibility','off');

txtslice = uicontrol('Style','text',...
        'Units','normalized','Position',[0.48 0.08 0.3 0.05],...
        'String',['Slice: 1/' num2str(size(V,3))],'FontSize',14);
hslider = uicontrol('Style','slider','Enable','off',...
        'Min', 1, 'Max', size(V,3), 'Value', 1, 'SliderStep', [1/(size(V,3)-1), 10/(size(V,3)-1)],...
        'Units','normalized','Position',[0.43 0.05 0.4 0.05],...
        'Callback', @sliceSlider);

set(himg, 'CData', max(V,[],3));

% param panel
hp = uipanel('Title','Parameters','FontSize',14,...
             'Position',[.05 .5 .2 .4]);

uicontrol('Parent',hp,'Style','text','HorizontalAlignment','Left',...
        'Units','normalized','Position',[0.05 0.85 0.9 0.1],...
        'String','Slices','FontWeight','Bold','FontSize',14);
uicontrol('Parent',hp,'Style','text',...
        'Units','normalized','Position',[0.05 0.75 0.2 0.1],...
        'TooltipString', 'Pick the starting depth',...
        'String','from','FontSize',14);
txtsfrom = uicontrol('Parent',hp,'Style','edit',...
        'Units','normalized','Position',[0.25 0.75 0.25 0.1],'HorizontalAlignment','Left',...
        'String','1','FontSize',14,...
        'Callback',@checkslice_from);
uicontrol('Parent',hp,'Style','text',...
        'Units','normalized','Position',[0.5 0.75 0.2 0.1],...
        'TooltipString', 'Pick the Last depth',...
        'String','to','FontSize',14);
txtsto = uicontrol('Parent',hp,'Style','edit',...
        'Units','normalized','Position',[0.7 0.75 0.25 0.1],'HorizontalAlignment','Left',...
        'String',num2str(size(V,3)),'FontSize',14,...
        'Callback',@checkslice_to);

uicontrol('Parent',hp,'Style','text','HorizontalAlignment','Left',...
        'Units','normalized','Position',[0.05 0.6 0.9 0.1],...
        'String','Image Intensity','FontWeight','Bold','FontSize',14);
uicontrol('Parent',hp,'Style','text',...
        'Units','normalized','Position',[0.05 0.5 0.2 0.1],...
        'TooltipString', 'Lowest input intensity so this or lower value will be zero in the output.',...
        'String','Low:','FontSize',14);
txtlow_in = uicontrol('Parent',hp,'Style','edit',...
        'Units','normalized','Position',[0.25 0.5 0.25 0.1],'HorizontalAlignment','Left',...
        'String','0.0','FontSize',14,...
        'Callback',@checklow_in);
uicontrol('Parent',hp,'Style','text',...
        'Units','normalized','Position',[0.5 0.5 0.2 0.1],...
        'TooltipString', 'Highest input intensity so this or higher value will be one in the output.',...
        'String','High:','FontSize',14);
txthigh_in = uicontrol('Parent',hp,'Style','edit',...
        'Units','normalized','Position',[0.7 0.5 0.25 0.1],'HorizontalAlignment','Left',...
        'String','1.0','FontSize',14,...
        'Callback',@checkhigh_in);

uicontrol('Parent',hp,'Style','pushbutton','String','Update','FontSize',14,...
        'Units','normalized','Position',[0.3 0.1 0.4 0.15],...
        'TooltipString', 'Update parameters and show the pre-processed image stack.',...
        'Callback',@updateCallback);
    
uicontrol('Style','pushbutton','String','Done','FontSize',14,...
        'Units','normalized','Position',[.1 .35 .1 .1],...
        'Callback',@traceCallback);

uiwait(gcf);
close(f);

    function checkslice_from(~,~)
        slice_fr = str2double(txtsfrom.String);
        slice_to = str2double(txtsto.String);
        if isnan(slice_fr) || slice_fr < 1 || slice_fr > slice_to,
            errordlg(['Starting slice must be an integer between 1 and ' txtsto.String '.']);
            txtsfrom.String = '1';
        end
    end

    function checkslice_to(~,~)
        slice_fr = str2double(txtsfrom.String);
        slice_to = str2double(txtsto.String);
        if isnan(slice_to) || slice_to > size(V,3) || slice_fr > slice_to,
            errordlg(['Last slice must be an integer between ' txtsfrom.String ' and ' num2str(size(V,3)) '.']);
            txtsto.String = num2str(size(V,3));
        end
    end

    function checklow_in(~,~)
        low_in = str2double(txtlow_in.String);
        high_in = str2double(txthigh_in.String);
        if isnan(low_in) || low_in < 0 || low_in > high_in,
            errordlg(['Lowest intensity must be a real number between 0 and ' txthigh_in.String '.']);
            txtlow_in.String = '0';
        end
    end

    function checkhigh_in(~,~)
        low_in = str2double(txtlow_in.String);
        high_in = str2double(txthigh_in.String);
        if isnan(high_in) || high_in > 1 || low_in > high_in,
            errordlg(['Highest intensity must be a real number between ' txtlow_in.String ' and 1.']);
            txthigh_in.String = '1';
        end
    end

    function viewSelectCallback(~,~)
        if h2d.Value
            hslider.Enable = 'off';
            set(himg, 'CData', max(Vproc,[],3));
        elseif h3d.Value
            hslider.Enable = 'on';
            set(himg, 'CData', Vproc(:,:,hslider.Value));
        end
    end

    function sliceSlider(~,~)
        hslider.Value = round(hslider.Value);
        set(himg, 'CData', Vproc(:,:,hslider.Value));
        txtslice.String = ['Slice: ' num2str(hslider.Value) '/' num2str(size(Vproc,3))];
    end

    function updateCallback(~,~)
        slice_fr = str2double(txtsfrom.String);
        slice_to = str2double(txtsto.String);
        low_in = str2double(txtlow_in.String);
        high_in = str2double(txthigh_in.String);
        
        fwait = waitbar(0.5, 'Pre-processing image stack');
        Vproc = im2double(V(:,:,slice_fr:slice_to));
        for z = 1:size(Vproc,3)
            Vproc(:,:,z) = imadjust(uint8(Vproc(:,:,z)),[low_in, high_in],[0,1]);
        end
        close(fwait);
        
        hslider.Max = size(Vproc,3);
        hslider.SliderStep = [1/(size(Vproc,3)-1), 10/size(Vproc,3)];
        hslider.Value = 1;
        txtslice.String = ['Slice: 1/' num2str(size(Vproc,3))];
        viewSelectCallback();
    end

    function traceCallback(~,~)
        slice_fr = str2double(txtsfrom.String);
        slice_to = str2double(txtsto.String);
        low_in = str2double(txtlow_in.String);
        high_in = str2double(txthigh_in.String);
        params = struct('sfr', slice_fr, ...
                'sto', slice_to, ...
                'low_in', low_in, ...
                'high_in', high_in);
        uiresume(gcbf);
    end
end

