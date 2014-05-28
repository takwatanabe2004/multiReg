function varargout = tak_gui_subplot_volumes(varargin)
%% varargout = tak_gui_subplot_volumes(varargin), the final entry must contain [] or slider spacing prooperty
%% varargin{1...} = insert as many 3d volume slices you want to display
%% varargin{end} = sliderspacing property
% 
% Example
% tak_gui_subplot_volumes(Iref,Ireg,warpgrid,positive_Jac,[min(Iref(:)) max(Iref(:))], [])
% 
% % axial
% tak_gui_subplot_volumes(Iref_fid2,Ireg_fid,Ihol_fid,Iref,[0,maxcmap],[1 fid.x(2)-fid.x(1)])
% 
% % saggital
% tak_gui_subplot_volumes(permute(Iref_fid2,[2 3 1]),...
%                     permute(Ireg_fid, [2 3 1]),...
%                     permute(Ihol_fid, [2 3 1]),...
%                     permute(Iref,[2 3 1]),...
%                     [0,maxcmap],[1 fid.x(2)-fid.x(1)])
% 
% % coronal
% tak_gui_subplot_volumes(permute(Iref_fid2,[1 3 2]),...
%                     permute(Ireg_fid, [1 3 2]),...
%                     permute(Ihol_fid, [1 3 2]),...
%                     permute(Iref,[1 3 2]),...
%                     [0,maxcmap],[1 fid.x(2)-fid.x(1)])  

% TAK_GUI_SUBPLOT_VOLUMES M-file for tak_gui_subplot_volumes.fig
%      TAK_GUI_SUBPLOT_VOLUMES, by itself, creates a new TAK_GUI_SUBPLOT_VOLUMES or raises the existing
%      singleton*.
%
%      H = TAK_GUI_SUBPLOT_VOLUMES returns the handle to a new TAK_GUI_SUBPLOT_VOLUMES or the handle to
%      the existing singleton*.
%
%      TAK_GUI_SUBPLOT_VOLUMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TAK_GUI_SUBPLOT_VOLUMES.M with the given input arguments.
%
%      TAK_GUI_SUBPLOT_VOLUMES('Property','Value',...) creates a new TAK_GUI_SUBPLOT_VOLUMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tak_gui_subplot_volumes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tak_gui_subplot_volumes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tak_gui_subplot_volumes

% Last Modified by GUIDE v2.5 06-May-2010 02:40:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tak_gui_subplot_volumes_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_subplot_volumes_OutputFcn, ...
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


% --- Executes just before tak_gui_subplot_volumes is made visible.
function tak_gui_subplot_volumes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tak_gui_subplot_volumes (see VARARGIN)

% Choose default command line output for tak_gui_subplot_volumes
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tak_gui_subplot_volumes wait for user response (see UIRESUME)
% uiwait(handles.figure1);

data.fill = 0;

data.size = size(varargin{1});

if ~isempty(varargin{size(varargin,2)})
    data.spacing1_small = varargin{size(varargin,2)}(1)/data.size(3);
    data.spacing1_big   = varargin{size(varargin,2)}(2)/data.size(3);   
    % Set slider value range accordingly
    set(handles.slider1,'SliderStep',[data.spacing1_small data.spacing1_big])
end

data.naxis = size(varargin,2)-1;

data.subplot_nrow = 1;  
if data.naxis >= 4    
    data.subplot_nrow = 2;    
end
data.subplot_ncol = ceil(data.naxis/data.subplot_nrow);

for idx = 1:data.naxis
    subplot(data.subplot_nrow,data.subplot_ncol,idx)
    data.images{idx} = varargin{idx};
    imagesc(data.images{idx}(:,:,1)'),axis('xy', 'image','off')
    if data.fill == 1
        axis('normal')
    end
end
colormap(gray)
setMyData(data)

% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_subplot_volumes_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data = getMyData();
sliderpos = get(hObject,'Value');
pos = round(data.size(3)*sliderpos);
if pos == 0; pos = 1; end
% keyboard
    for idx = 1:data.naxis
        subplot(data.subplot_nrow,data.subplot_ncol,idx)
        imagesc(data.images{idx}(:,:,pos)'),axis('xy', 'image','off')
        if data.fill == 1
            axis('normal')
        end
        if idx ==1
            title(strcat(num2str(pos),{' out ot '},num2str(data.size(3))))
        end
    end
setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function data=getMyData()
data = getappdata(gca,'showdata');

function setMyData(data)
setappdata(gca,'showdata',data);
