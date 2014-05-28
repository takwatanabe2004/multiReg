function varargout = tak_gui_sweep_onevar_image(varargin)
%% varargin{1} = the image
%% varargin{2} = the function handle
%% varargin{3} = the sweep range of the variable of interest
% TAK_GUI_SWEEP_ONEVAR_IMAGE M-file for tak_gui_sweep_onevar_image.fig
%      TAK_GUI_SWEEP_ONEVAR_IMAGE, by itself, creates a new TAK_GUI_SWEEP_ONEVAR_IMAGE or raises the existing
%      singleton*.
%
%      H = TAK_GUI_SWEEP_ONEVAR_IMAGE returns the handle to a new TAK_GUI_SWEEP_ONEVAR_IMAGE or the handle to
%      the existing singleton*.
%
%      TAK_GUI_SWEEP_ONEVAR_IMAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TAK_GUI_SWEEP_ONEVAR_IMAGE.M with the given input arguments.
%
%      TAK_GUI_SWEEP_ONEVAR_IMAGE('Property','Value',...) creates a new TAK_GUI_SWEEP_ONEVAR_IMAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tak_gui_sweep_onevar_image_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tak_gui_sweep_onevar_image_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tak_gui_sweep_onevar_image

% Last Modified by GUIDE v2.5 05-May-2010 02:01:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tak_gui_sweep_onevar_image_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_sweep_onevar_image_OutputFcn, ...
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


% --- Executes just before tak_gui_sweep_onevar_image is made visible.
function tak_gui_sweep_onevar_image_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tak_gui_sweep_onevar_image (see VARARGIN)

% Choose default command line output for tak_gui_sweep_onevar_image
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% load input into structure 'data'
data.image = varargin{1};
data.sizes = size(data.image);

data.fhandle = varargin{2};

data.var = varargin{3}; % input should be sweep range of variable
data.var_low = data.var(1);
data.var_hi = data.var(end);
% the slider spacing increment must be between 0 and 1...so it is proportionate
data.var_spacing1 = (data.var(2) - data.var(1))/(data.var_hi - data.var_low);

% Set slider value range accordingly
set(handles.slider1,'Max',data.var_hi,'Min',data.var_low,'SliderStep',[data.var_spacing1 0.05]) 

% For default, use the halfway point of sweeprange
data.slider1_pos = 0.5 * (data.var_hi - data.var_low) + data.var_low; 
set(handles.slider1,'value',data.slider1_pos); % set slider location

% Display result
set(gcf,'CurrentAxes',handles.axes1)
imagesc_tak2(data.fhandle(data.slider1_pos)')
title(strcat({'Variable = '},num2str(data.slider1_pos),{' (Range: '},num2str(data.var_low),{' to '},num2str(data.var_hi),')'))
colorbar('location','northoutside')
% update data (also necessary to be able to get content of 'data' in
% function call from slider_callback
setMyData(data)
% UIWAIT makes tak_gui_sweep_onevar_image wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_sweep_onevar_image_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Retrieve data
data = getMyData();
data.slider1_pos = get(hObject,'Value');

% display result
set(gcf,'CurrentAxes',handles.axes1)
imagesc_tak2(data.fhandle(data.slider1_pos)')
title(strcat({'Variable = '},num2str(data.slider1_pos),{' (Range: '},num2str(data.var_low),{' to '},num2str(data.var_hi),')'))
colorbar('location','northoutside')
% update data
setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function data=getMyData()
data = getappdata(gca,'showdata'); 

function setMyData(data)
setappdata(gca,'showdata',data); 
