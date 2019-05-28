function varargout = view3Dmodel(varargin)
% VIEW3DMODEL MATLAB code for view3Dmodel.fig
%      VIEW3DMODEL, by itself, creates a new VIEW3DMODEL or raises the existing
%      singleton*.
%
%      H = VIEW3DMODEL returns the handle to a new VIEW3DMODEL or the handle to
%      the existing singleton*.
%
%      VIEW3DMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW3DMODEL.M with the given input arguments.
%
%      VIEW3DMODEL('Property','Value',...) creates a new VIEW3DMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view3Dmodel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view3Dmodel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view3Dmodel

% Last Modified by GUIDE v2.5 15-Dec-2017 10:59:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view3Dmodel_OpeningFcn, ...
                   'gui_OutputFcn',  @view3Dmodel_OutputFcn, ...
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


% --- Executes just before view3Dmodel is made visible.
function view3Dmodel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view3Dmodel (see VARARGIN)

% Choose default command line output for view3Dmodel
handles.output = hObject;

handles.dx = 4;

col1 = [235 90  235]/255; % ddaD color
col2 = [80  200 80 ]/255; % ddaE color

handles.V = varargin{1};
handles.track_time = varargin{2};
handles.swcs = varargin{3};
handles.soma_pos = varargin{4};
handles.configs = varargin{5};
handles.spline = varargin{6};
handles.zflip = varargin{7};

[H, W, D, T] = size(handles.V);
x = round(W/2); y = round(H/2); z = round(D/2); t = round(T/2);
set(handles.xslider, 'Min', 1, 'Max', W, 'Value', x, 'SliderStep', [1/(W-1), 10/(W-1)]);
set(handles.xval, 'String', x);
set(handles.yslider, 'Min', 1, 'Max', H, 'Value', y, 'SliderStep', [1/(H-1), 10/(H-1)]);
set(handles.yval, 'String', t);
set(handles.zslider, 'Min', 1, 'Max', D, 'Value', z, 'SliderStep', [1/(D-1), 2/(D-1)]);
set(handles.zval, 'String', z);
set(handles.tslider, 'Min', 1, 'Max', T, 'Value', 1, 'SliderStep', [1/(T-1), 10/(T-1)]);
set(handles.tval, 'String', handles.track_time);
set(handles.xmax, 'String', W);
set(handles.ymax, 'String', H);
set(handles.zmax, 'String', D);
set(handles.tmax, 'String', T+handles.track_time-1);

% draw image stack as three slices
axes(handles.axes1)
[xx, yy, zz] = meshgrid(x, 1:H, 1:D);
cc = handles.V(:,x,:,t);
handles.hx = surface(squeeze(xx), squeeze(yy), squeeze(zz), squeeze(cc),...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');
handles.hx.AmbientStrength = 1;
handles.hx.SpecularStrength = 0;

[xx, yy, zz] = meshgrid(1:W, y, 1:D);
cc = handles.V(y,:,:,t);
handles.hy = surface(squeeze(xx), squeeze(yy), squeeze(zz), squeeze(cc),...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');
handles.hy.AmbientStrength = 1;
handles.hy.SpecularStrength = 0;

[xx, yy, zz] = meshgrid(1:W, 1:H, z);
cc = handles.V(:,:,z,t);
handles.hz = surface(xx, yy, zz, cc,...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');
handles.hz.AmbientStrength = 1;
handles.hz.SpecularStrength = 0;

colormap(gray(256));
view(3)
% axis equal tight



% locate soma
ddad_somaX = round(handles.swcs{1}(1,3));
ddad_somaY = round(handles.swcs{1}(1,4));
ddae_somaX = round(handles.swcs{2}(1,3));
ddae_somaY = round(handles.swcs{2}(1,4));
somaX = round((ddad_somaX + ddae_somaX)/2);
somaY = round((ddad_somaY + ddae_somaY)/2);

% create baseline mask
sizeI = [size(handles.V,1), size(handles.V,2)];
I_mask1 = swc2pixel( handles.swcs{1}, sizeI );
I_mask2 = swc2pixel( handles.swcs{2}, sizeI );
I2D_mask = imdilate(I_mask1 | I_mask2, ones(3));

% get bounding box
[ii, jj] = find(I2D_mask);
xrange = [min(jj) max(jj)];
yrange = [min(ii) max(ii)];
handles.I_mask1_crop = I_mask1(yrange(1):yrange(2),xrange(1):xrange(2));
handles.I_mask2_crop = I_mask2(yrange(1):yrange(2),xrange(1):xrange(2));

% init FFD
[Ox_bs, Oz_bs] = ffd_init_from_config(handles.configs(handles.track_time+1,:), ...
        somaX, handles.dx, handles.zflip{handles.track_time+1});

% transform dendrite
Tx = ffd_interpolate(Ox_bs, handles.spline);
Tz = ffd_interpolate(Oz_bs, handles.spline)+1;
[cx, cy] = meshgrid(round(Tx), yrange(1):yrange(2));
[cz, ~] = meshgrid(round(Tz), yrange(1):yrange(2));
ddad_somaZ = cz(ddad_somaY - yrange(1) + 1, ddad_somaX - xrange(1) + 1);
ddae_somaZ = cz(ddae_somaY - yrange(1) + 1, ddae_somaX - xrange(1) + 1);
cx = cx(:); cy = cy(:); cz = cz(:);

% create 3D surface of ddaD
sizeI3D = [sizeI 1+range(cz)];
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(handles.I_mask1_crop), cx(handles.I_mask1_crop), cz(handles.I_mask1_crop))) = 1;
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
[xx3, yy3, zz3] = meshgrid(1:size(I3D,2), 1:size(I3D,1), 1:size(I3D,3));
I3D((xx3 - ddad_somaX - 2).^2 + (yy3 - ddad_somaY - 2).^2 + ...
        (zz3 - ddad_somaZ - 7).^2 < handles.swcs{1}(1,6)^2) = 1;
sf1 = isosurface(I3D, .5);
% relocate soma to (0,0,1)
sf1.vertices = bsxfun(@minus, sf1.vertices, [2, 2, 7 + ddad_somaZ - 6]);
sf1.vertices(:,3) = sf1.vertices(:,3) * 0.624/1.06;

% create 3D surface of ddaE
I3D = zeros(sizeI3D);
I3D(sub2ind(sizeI3D, cy(handles.I_mask2_crop), cx(handles.I_mask2_crop), cz(handles.I_mask2_crop))) = 1;
I3D = padarray(I3D, [2 2 7], 0, 'both');
I3D = imdilate(I3D, ones(3,3,3));
I3D((xx3 - ddae_somaX - 2).^2 + (yy3 - ddae_somaY - 2).^2 + ...
        (zz3 - ddae_somaZ - 7).^2 < handles.swcs{2}(1,6)^2) = 1;
sf2 = isosurface(I3D, .5);
% relocate soma to (0,0,1)
sf2.vertices = bsxfun(@minus, sf2.vertices, [2, 2, 7 + ddae_somaZ - 6]);
sf2.vertices(:,3) = sf2.vertices(:,3) * 0.624/1.06;

% draw the neuron 3D model as surfaces
handles.hmodel1 = patch(sf1);
handles.hmodel1.FaceColor = col1;
handles.hmodel1.EdgeColor = 'none';

handles.hmodel2 = patch(sf2);
handles.hmodel2.FaceColor = col2;
handles.hmodel2.EdgeColor = 'none';

daspect([1 1 .3])
camproj('perspective')
camlight 
lighting gouraud
set(handles.axes1,'Zdir','reverse');

handles.somaY = somaY;
handles.ddad_somaX = ddad_somaX;
handles.ddad_somaY = ddad_somaY;
handles.ddae_somaX = ddae_somaX;
handles.ddae_somaY = ddae_somaY;
handles.xrange = xrange;
handles.yrange = yrange;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes view3Dmodel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view3Dmodel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function xslider_Callback(hObject, eventdata, handles)
% hObject    handle to xslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
x = round(get(hObject,'Value'));
set(handles.xval, 'String', x);
t = round(get(handles.tslider, 'Value'));
updateX(x, t, handles.hx, handles.V);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function xslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function yslider_Callback(hObject, eventdata, handles)
% hObject    handle to yslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
y = round(get(hObject,'Value'));
set(handles.yval, 'String', y);
t = round(get(handles.tslider, 'Value'));
updateY(y, t, handles.hy, handles.V);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zslider_Callback(hObject, eventdata, handles)
% hObject    handle to zslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
z = round(get(hObject,'Value'));
set(handles.zval, 'String', z);
t = round(get(handles.tslider, 'Value'));
updateZ(z, t, handles.hz, handles.V);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function zslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function tslider_Callback(hObject, eventdata, handles)
% hObject    handle to tslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
t = round(get(hObject,'Value'));
set(handles.tval, 'String', t+handles.track_time-1);
% update slices
x = round(get(handles.xslider, 'Value'));
y = round(get(handles.yslider, 'Value'));
z = round(get(handles.zslider, 'Value'));
updateX(x, t, handles.hx, handles.V);
updateY(y, t, handles.hy, handles.V);
updateZ(z, t, handles.hz, handles.V);
% update 3D model


if get(handles.neuronvis, 'Value')
% transform dendrite in Z-axis
    sizeI = [size(handles.V,1) size(handles.V,2)];
    t = t + handles.track_time;
    ddadX = handles.soma_pos{1}(t,1);
    ddadY = handles.soma_pos{1}(t,2);
    ddaeX = handles.soma_pos{2}(t,1);
    ddaeY = handles.soma_pos{2}(t,2);
    newsomaX = (ddadX + ddaeX)/2;
    dsomaY = round(handles.somaY - (ddadY+ddaeY)/2);

    [Ox,Oz] = ffd_init_from_config(handles.configs(t,:), newsomaX, handles.dx, handles.zflip{t});
    
    % transform dendrite in 2D
    Tx = ffd_interpolate(Ox, handles.spline);
    newyrange = handles.yrange - dsomaY;
    [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
    cx = cx(:); cy = cy(:);
    idxIn = cy >= 1 & cy <= sizeI(1) & cx >= 1 & cx <= sizeI(2);

    Tz = ffd_interpolate(Oz, handles.spline)+1;
    [cz, ~] = meshgrid(round(Tz) - min(round(Tz))+1, newyrange(1):newyrange(2));
    ddadZ = cz(handles.ddad_somaY - handles.yrange(1) + 1, handles.ddad_somaX - handles.xrange(1) + 1);
    ddaeZ = cz(handles.ddae_somaY - handles.yrange(1) + 1, handles.ddae_somaX - handles.xrange(1) + 1);
    cz = cz(:);
    sizeI3D = [sizeI 1+range(cz)];
    % update 3D dendrite ddaD
    I3D = zeros(sizeI3D);
    I3D(sub2ind(sizeI3D, cy(handles.I_mask1_crop(:)&idxIn), ...
            cx(handles.I_mask1_crop(:)&idxIn), cz(handles.I_mask1_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    [xx3, yy3, zz3] = meshgrid(1:size(I3D,2), 1:size(I3D,1), 1:size(I3D,3));
    I3D((xx3 - ddadX - 2).^2 + (yy3 - ddadY - 2).^2 + ...
            (zz3 - ddadZ - 7).^2 < handles.swcs{1}(1,6)^2) = 1;
    sf1 = isosurface(I3D, .5);
    % relocate soma to (0,0,1)
    handles.hmodel1.Vertices = bsxfun(@minus, sf1.vertices, [2, 2, 7 + ddadZ - 6]);
    handles.hmodel1.Vertices(:,3) = handles.hmodel1.Vertices(:,3) * 0.624/1.06;
    handles.hmodel1.Faces = sf1.faces;
    % update 3D dendrite ddaE
    I3D = zeros(sizeI3D);
    I3D(sub2ind(sizeI3D, cy(handles.I_mask2_crop(:)&idxIn), ...
            cx(handles.I_mask2_crop(:)&idxIn), cz(handles.I_mask2_crop(:)&idxIn))) = 1;
    I3D = padarray(I3D, [2 2 7], 0, 'both');
    I3D = imdilate(I3D, ones(3,3,3));
    I3D((xx3 - ddaeX - 2).^2 + (yy3 - ddaeY - 2).^2 + ...
            (zz3 - ddaeZ - 7).^2 < handles.swcs{2}(1,6)^2) = 1;
    
    sf2 = isosurface(I3D, .5);
    % relocate soma to (0,0,1)
    handles.hmodel2.Vertices = bsxfun(@minus, sf2.vertices, [2, 2, 7 + ddaeZ - 6]);
    handles.hmodel2.Vertices(:,3) = handles.hmodel2.Vertices(:,3) * 0.624/1.06;
    handles.hmodel2.Faces = sf2.faces;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in neuronvis.
function neuronvis_Callback(hObject, eventdata, handles)
% hObject    handle to neuronvis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuronvis
vis = get(hObject, 'Value');
if vis
    handles.hmodel1.Visible = 'on';
    handles.hmodel2.Visible = 'on';
else
    handles.hmodel1.Visible = 'off';
    handles.hmodel2.Visible = 'off';
end
% Update handles structure
guidata(hObject, handles);

% --- auxiliary functions
function updateX(x, t, hx, V)
    cc = V(:,x,:,t);
    set(hx, 'CData', squeeze(cc));
    xx = get(hx, 'XData');
    xx(:) = x;
    set(hx, 'XData', xx);
    
function updateY(y, t, hy, V)
    cc = V(y,:,:,t);
    set(hy, 'CData', squeeze(cc));
    yy = get(hy, 'YData');
    yy(:) = y;
    set(hy, 'YData', yy);
    
function updateZ(z, t, hz, V)
    cc = V(:,:,z,t);
    set(hz, 'CData', cc);
    zz = get(hz, 'ZData');
    zz(:) = z;
    set(hz, 'ZData', zz);
