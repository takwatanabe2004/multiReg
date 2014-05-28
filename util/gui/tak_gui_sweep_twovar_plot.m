function varargout = tak_gui_sweep_twovar_plot(varargin)
%% varargin{1} = the independent variable
%% varargin{2} = the function handle
%% varargin{3} = the sweep range of the 1st variable of interest
%% varargin{4} = the sweep range of the 1st variable of interest
%% varargin{5} = optional dependent variable you want to subplot against
%% (in the future will allow any amount of dep var...not too hard to do)
% TAK_GUI_SWEEP_TWOVAR_PLOT M-file for tak_gui_sweep_twovar_plot.fig
%      TAK_GUI_SWEEP_TWOVAR_PLOT, by itself, creates a new TAK_GUI_SWEEP_TWOVAR_PLOT or raises the existing
%      singleton*.
%
%      H = TAK_GUI_SWEEP_TWOVAR_PLOT returns the handle to a new TAK_GUI_SWEEP_TWOVAR_PLOT or the handle to
%      the existing singleton*.
%
%      TAK_GUI_SWEEP_TWOVAR_PLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TAK_GUI_SWEEP_TWOVAR_PLOT.M with the given input arguments.
%
%      TAK_GUI_SWEEP_TWOVAR_PLOT('Property','Value',...) creates a new TAK_GUI_SWEEP_TWOVAR_PLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tak_gui_sweep_twovar_plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tak_gui_sweep_twovar_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tak_gui_sweep_twovar_plot

% Last Modified by GUIDE v2.5 05-May-2010 02:39:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tak_gui_sweep_twovar_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_sweep_twovar_plot_OutputFcn, ...
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


% --- Executes just before tak_gui_sweep_twovar_plot is made visible.
function tak_gui_sweep_twovar_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tak_gui_sweep_twovar_plot (see VARARGIN)

% Choose default command line output for tak_gui_sweep_twovar_plot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% load input into structure 'data'
data.x = varargin{1}; % the independent variable (x-axis)
data.fhandle = varargin{2};
data.var1 = varargin{3}; % input should be sweep range of variable
data.var2 = varargin{4}; % input should be sweep range of variable
if(size(varargin,2)==5)
    data.y1 = varargin{5};
else
    data.y1 = [];
end

% Setup 1st variable wrt slider1
data.var1_low = data.var1(1);
data.var1_hi = data.var1(end);
% the slider spacing increment must be between 0 and 1...so it is proportionate
data.var1_spacing = (data.var1(2) - data.var1(1))/(data.var1_hi - data.var1_low);
% Set slider value range accordingly
set(handles.slider1,'Max',data.var1_hi,'Min',data.var1_low,'SliderStep',[data.var1_spacing 0.05]) 
% For default, use the halfway point of sweeprange
data.slider1_pos = 0.5 * (data.var1_hi - data.var1_low) + data.var1_low; 
set(handles.slider1,'value',data.slider1_pos); % set slider location

% Setup 2nd variable wrt slider2
data.var2_low = data.var2(1);
data.var2_hi = data.var2(end);
% the slider spacing increment must be between 0 and 1...so it is proportionate
data.var2_spacing = (data.var2(2) - data.var2(1))/(data.var2_hi - data.var2_low);
% Set slider value range accordingly
set(handles.slider2,'Max',data.var2_hi,'Min',data.var2_low,'SliderStep',[data.var2_spacing 0.05]) 
% For default, use the halfway point of sweeprange
data.slider2_pos = 0.5 * (data.var2_hi - data.var2_low) + data.var2_low; 
set(handles.slider2,'value',data.slider2_pos); % set slider location

% Display result
plot(handles.axes1,data.x,data.fhandle(data.slider1_pos,data.slider2_pos))
if(~isempty(data.y1)) hold on; plot(data.x,data.y1,'r'); end
% set(handles.axes1,'XMinorTick','on')
grid on
part1 = strcat({'Var1 = '},num2str(data.slider1_pos),{' ('},num2str(data.var1_low),{' to '},num2str(data.var1_hi),{'), '});
part2 = strcat({'Var2 = '},num2str(data.slider2_pos),{' ('},num2str(data.var2_low),{' to '},num2str(data.var2_hi),')');
title(strcat(part1,part2))
% update data (also necessary to be able to get content of 'data' in
% function call from slider_callback
setMyData(data)
% UIWAIT makes tak_gui_sweep_twovar_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_sweep_twovar_plot_OutputFcn(hObject, eventdata, handles) 
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

% Display result
hold off
plot(handles.axes1,data.x,data.fhandle(data.slider1_pos,data.slider2_pos))
if(~isempty(data.y1)) hold on; plot(data.x,data.y1,'r'); end
% set(handles.axes1,'XMinorTick','on')
grid on
part1 = strcat({'Var1 = '},num2str(data.slider1_pos),{' ('},num2str(data.var1_low),{' to '},num2str(data.var1_hi),{'), '});
part2 = strcat({'Var2 = '},num2str(data.slider2_pos),{' ('},num2str(data.var2_low),{' to '},num2str(data.var2_hi),')');
title(strcat(part1,part2))
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

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Retrieve data
data = getMyData();
data.slider2_pos = get(hObject,'Value');

% Display result
hold off
plot(handles.axes1,data.x,data.fhandle(data.slider1_pos,data.slider2_pos))
if(~isempty(data.y1)) hold on; plot(data.x,data.y1,'r'); end
% set(handles.axes1,'XMinorTick','on')
grid on
part1 = strcat({'Var1 = '},num2str(data.slider1_pos),{' ('},num2str(data.var1_low),{' to '},num2str(data.var1_hi),{'), '});
part2 = strcat({'Var2 = '},num2str(data.slider2_pos),{' ('},num2str(data.var2_low),{' to '},num2str(data.var2_hi),')');
title(strcat(part1,part2))
% update data
setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
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
