function varargout = GUI_edit_trace(varargin)
% GUI_EDIT_TRACE MATLAB code for GUI_edit_trace.fig
%      GUI_EDIT_TRACE, by itself, creates a new GUI_EDIT_TRACE or raises the existing
%      singleton*.
%
%      H = GUI_EDIT_TRACE returns the handle to a new GUI_EDIT_TRACE or the handle to
%      the existing singleton*.
%
%      GUI_EDIT_TRACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_EDIT_TRACE.M with the given input arguments.
%
%      GUI_EDIT_TRACE('Property','Value',...) creates a new GUI_EDIT_TRACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_edit_trace_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_edit_trace_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_edit_trace

% Last Modified by GUIDE v2.5 15-Dec-2016 12:17:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_edit_trace_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_edit_trace_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_edit_trace is made visible.
function GUI_edit_trace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_edit_trace (see VARARGIN)

% Choose default command line output for GUI_edit_trace
handles.output = hObject;

addpath(genpath('matlab_bgl'));

% ask user to select image stack and swc file to edit
% load('matlab.mat');
img_dir = uigetdir('~/Desktop/LipingsData/SEM_ Tiled C4da neurons/', 'Select Image Stack Folder to Open');
if img_dir == 0
    errordlg('No Image Stack Directory Selected');
    delete(handles.figure1)
    return;
end
[swcfilename,swcpathname,~] = uigetfile('*.swc', 'Select SWC File to Open');
if swcfilename == 0
    errordlg('No SWC File Selected');
    delete(handles.figure1)
    return;
end

% get max intensity from user (default value is 150)
prompt = {'Enter Max Intensity (0-255):'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'150'};
dinp = inputdlg(prompt,dlg_title,num_lines,defaultans);
max_intensity = min(max(str2double(dinp{1}), 0), 255);

handles.V = IV_loadimg([img_dir filesep]);
handles.swcdata = read_swc_file([swcpathname swcfilename]);
% save('matlab.mat', 'img_dir', 'swcpathname', 'swcfilename');

% undo option
handles.ori_swcdata = handles.swcdata(:,7);
handles.action = zeros(size(handles.ori_swcdata));
handles.action_num = 0;

% update slider and axis
axes(handles.img_axes)
handles.img_hdl = imshow(handles.V(:,:,1), [0 max_intensity]);
set(handles.img_hdl,'HitTest','off');

hold on
handles.dendrite_plot = [];
nb = handles.swcdata(:,[1 7]);
nb(any(nb < 1,2),:) = [];
G = sparse(nb(:,1), nb(:,2), 1, size(handles.swcdata,1), size(handles.swcdata,1));
G = G + G';
[cc, sizes] = components(G);
col = jet(max(5, length(sizes)+1));
for j = 1:length(sizes)
    idx = all(cc(nb)==j,2);
    xx = [handles.swcdata(nb(idx,1),3), handles.swcdata(nb(idx,2),3), nan(sum(idx),1)]';
    yy = [handles.swcdata(nb(idx,1),4), handles.swcdata(nb(idx,2),4), nan(sum(idx),1)]';
    handles.dendrite_plot(j) = plot(xx(:), yy(:), 'LineWidth', 1, 'Color', col(j,:));
end
hold off

FNUM = size(handles.V,3);
set(handles.depth_label, 'String', [num2str(0) ' / ' num2str(FNUM)]);
set(handles.depth_slider, 'Min', 1, 'Max', FNUM, 'Value', 1, 'SliderStep', [1/(FNUM-1), 10/(FNUM-1)]);

handles.button_down = false;
handles.cutline = zeros(2,2);
handles.cutline_plot = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_edit_trace wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_edit_trace_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    varargout{1} = handles.output;
catch
end


% --- Executes on slider movement.
function depth_slider_Callback(hObject, eventdata, handles)
% hObject    handle to depth_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
cur_depth = round(get(hObject,'Value'));
set(hObject,'Value',cur_depth);
FNUM = size(handles.V,3);
set(handles.img_hdl, 'CData', handles.V(:,:,cur_depth));
set(handles.depth_label, 'String', [num2str(cur_depth-1) ' / ' num2str(FNUM)]);


% --- Executes during object creation, after setting all properties.
function depth_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in undo_button.
function undo_button_Callback(hObject, eventdata, handles)
% hObject    handle to undo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% undo splitting
action_idx = handles.action == handles.action_num;
handles.swcdata(action_idx,7) = handles.ori_swcdata(action_idx);
handles.action(action_idx) = 0;

% update trace
axes(handles.img_axes)
handles.dendrite_plot = redraw_trace( handles.swcdata, handles.dendrite_plot );

% update order of action
handles.action_num = handles.action_num - 1;
if handles.action_num == 0
    set(handles.undo_button, 'Enable', 'off');
end

guidata(hObject, handles);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[savename, savepath] = uiputfile('*.swc', 'Save Trace File');
if savename ~= 0
    savefile = [savepath, savename];
    fileID = fopen(savefile, 'w');
    fprintf(fileID,'%d %d %.3f %.3f %.3f %.4f %d\n',handles.swcdata');
    fclose(fileID);
    msgbox(['Save SWC file: ' savefile ' completed'], 'Save Success');
else
    msgbox('Save SWC cancelled');
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.button_down
    cursorPoint = get(handles.img_axes, 'CurrentPoint');
    curX = cursorPoint(1,1);
    curY = cursorPoint(1,2);

    xLimits = get(handles.img_axes, 'xlim');
    yLimits = get(handles.img_axes, 'ylim');

    if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits))
        %disp(['Cursor coordinates are (' num2str(curX) ', ' num2str(curY) ').']);
        handles.button_down = true;
        
        axes(handles.img_axes)
        hold on
        handles.cutline = repmat([curX, curY], 2, 1);
        handles.cutline_plot = plot(handles.cutline(:,1), handles.cutline(:,2), 'r-o');
        hold off

        guidata(hObject, handles);
    else
        %disp('Cursor is outside bounds of image.');
    end
end


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.button_down
    cursorPoint = get(handles.img_axes, 'CurrentPoint');
    curX = cursorPoint(1,1);
    curY = cursorPoint(1,2);

    xLimits = get(handles.img_axes, 'xlim');
    yLimits = get(handles.img_axes, 'ylim');

    if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits))
        %disp(['Current Cursor coordinates are (' num2str(curX) ', ' num2str(curY) ').']);
        
        handles.cutline(2,1) = curX;
        handles.cutline(2,2) = curY;
        
        set(handles.cutline_plot, 'XData', handles.cutline(:,1));
        set(handles.cutline_plot, 'YData', handles.cutline(:,2));
        
        guidata(hObject, handles);
    else
        %disp('Current Cursor is outside bounds of image.');
    end
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.button_down
    %disp('Release Mouse')
    handles.button_down = false;
    
    [ cutindex ] = check_intersection( handles.swcdata, handles.cutline );
    
    % remove cutline from GUI
    delete(handles.cutline_plot);
    
    % keep track order of action
    action_idx = cutindex & handles.action == 0;
    if ~any(action_idx)
        guidata(hObject, handles);
        return;
    end
    if handles.action_num == 0
        set(handles.undo_button, 'Enable', 'on');
    end
    handles.action_num = handles.action_num + 1;
    %disp(handles.action_num)
    handles.action(action_idx) = handles.action_num;
    
    % split trace by cutline
    handles.swcdata(action_idx,7) = -1;
    
    % redraw trace plot
    axes(handles.img_axes)
    handles.dendrite_plot = redraw_trace( handles.swcdata, handles.dendrite_plot );
    
    guidata(hObject, handles);
end
