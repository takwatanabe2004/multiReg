function varargout = tak_gui_show_slices_cor_overlay(varargin)
%% This function displays the slices of a 3d volume from the 3-coordinates...i.e. the 
%  axial (xy), coronal (xz), and saggital (yz) view.
%  varargin{1} = image volume to display slices
%  varargin{2} = image volume to overlay
%  varargin{3} = overlay color transparency strength (default 0.25)
%  varargin{4} = clim = [lo hi] (defautl [])
%  varargin{5} = slider spacing setting (#s of index to slide) = [smallx bigx smally bigy smallz bigz].
%  varargin{6} = structure data (data.cbar, data.cmap, data.fill)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tak_gui_show_slices_cor_overlay_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_show_slices_cor_overlay_OutputFcn, ...
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

% --- Executes just before tak_gui_show_slices_cor_overlay is made visible.
function tak_gui_show_slices_cor_overlay_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tak_gui_show_slices_cor_overlay
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%---------- Process input data first ----------------------%
% main volume to display
data.voxelvolume = varargin{1};
data.sizes = size(data.voxelvolume);

% volume to overlay
data.voxelvolume2 = varargin{2};

data.color_strength = 0.25;
if (size(varargin,2)>=3) && (~isempty(varargin{3}))
    data.color_strength = varargin{3};
end

data.clims = [];    
if(size(varargin,2)>=4) && (~isempty(varargin{4}))
    data.clims = varargin{4};    
end

if(size(varargin,2)>=5) && (~isempty(varargin{5}))
    % Set slider value range accordingly
    % The order is in a hokie way since i wanted the input to be in the form  
    % [smallx bigx smally bigy smallz bigz]
    data.spacing2_small = varargin{5}(1)/data.sizes(2);
    data.spacing2_big   = varargin{5}(2)/data.sizes(2);
    set(handles.slider2,'SliderStep',[data.spacing2_small data.spacing2_big])
else
    data.spacing2_small = 1/data.sizes(2);
    data.spacing2_big   = 10/data.sizes(2);
    set(handles.slider2,'SliderStep',[data.spacing2_small data.spacing2_big])
end

data.cbar = 0; % if cbara ==1, display colorbar on top
data.cmap = 0; % if cmap == 1, show image with color
data.fill = 1; % if fill == 1, fill axis (ie axis('normal'))
if(size(varargin,2)==6) && (~isempty(varargin{6}))
    data.cbar = varargin{6}.cbar; 
    data.cmap = varargin{6}.cmap; 
    data.fill = varargin{6}.fill; 
end
%\-------------------------------------------------------------------------%

data.HandleWindow = gcf;
data.posy = 0.5; 

set(handles.slider2,'value',data.posy);

posy=round(data.sizes(2)*data.posy);

% Display Coronal view
set(gcf,'CurrentAxes',handles.axes2)
if(~isempty(data.clims))
    imagesc(squeeze(data.voxelvolume(:,posy,:))',[data.clims(1) data.clims(2)])
else    
    imagesc(squeeze(data.voxelvolume(:,posy,:))')
end
%-----------Include this part for the transparent overlay---------%
hold on
transparent = scale_uint8(  squeeze(data.voxelvolume2(:,posy,:))  )';  

% set transparanecy color in rgb.
colour = cat(3, 1*ones(size(transparent)), 0*ones(size(transparent)), 1*ones(size(transparent)));

data.hcolor_coronal = imagesc(colour); 
hold off
set(data.hcolor_coronal, 'AlphaData', transparent*data.color_strength)
%\----------------------------------------------------------------%
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Coronal View: yslice = '},num2str(posy), {' out of '} ,num2str(data.sizes(2))))
setMyData(data)

% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_show_slices_cor_overlay_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data = getMyData();
data.posy = get(handles.slider2,'Value');
posy=round(data.sizes(2)*data.posy);
if posy == 0; posy = 1; end

% Display Image Now
set(gcf,'CurrentAxes',handles.axes2)
if(~isempty(data.clims))
    imagesc(squeeze(data.voxelvolume(:,posy,:))',[data.clims(1) data.clims(2)])
else    
    imagesc(squeeze(data.voxelvolume(:,posy,:))')
end
%-----------Include this part for the transparent overlay---------%
hold on
transparent = scale_uint8(  squeeze(data.voxelvolume2(:,posy,:))  )';  

% set transparanecy color in rgb.
colour = cat(3, 1*ones(size(transparent)), 0*ones(size(transparent)), 1*ones(size(transparent)));

data.hcolor_coronal = imagesc(colour); 
hold off
set(data.hcolor_coronal, 'AlphaData', transparent*data.color_strength)
%\----------------------------------------------------------------%
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Coronal View: yslice = '},num2str(posy), {' out of '} ,num2str(data.sizes(2))))
setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function data=getMyData()
data=getappdata(gcf,'showdata');

function setMyData(data)
setappdata(data.HandleWindow,'showdata',data);
