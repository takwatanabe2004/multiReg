function varargout = tak_gui_show_slices_sag_overlay(varargin)
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
                   'gui_OpeningFcn', @tak_gui_show_slices_sag_overlay_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_show_slices_sag_overlay_OutputFcn, ...
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

% --- Executes just before tak_gui_show_slices_sag_overlay is made visible.
function tak_gui_show_slices_sag_overlay_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tak_gui_show_slices_sag_overlay
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
    data.spacing3_small = varargin{5}(1)/data.sizes(1);
    data.spacing3_big   = varargin{5}(2)/data.sizes(1);
    set(handles.slider3,'SliderStep',[data.spacing3_small data.spacing3_big])
else
    data.spacing3_small = 1/data.sizes(1);
    data.spacing3_big   = 10/data.sizes(1);
    set(handles.slider3,'SliderStep',[data.spacing3_small data.spacing3_big])
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
data.posx = 0.25; data.posy = 0.5; data.posz = 0.5;

% Now set the handles and figure according to the slider location
set(handles.slider3,'value',data.posz);

posx=round(data.sizes(1)*data.posx);

% Display Saggital view
set(gcf,'CurrentAxes',handles.axes3)
if(~isempty(data.clims))
    imagesc(squeeze(data.voxelvolume(posx,:,:))',[data.clims(1) data.clims(2)])
else
    imagesc(squeeze(data.voxelvolume(posx,:,:))')
end
%-----------Include this part for the transparent overlay---------%
hold on
transparent = scale_uint8(  squeeze(data.voxelvolume2(posx,:,:))  )';  

% set transparanecy color in rgb.
colour = cat(3, 1*ones(size(transparent)), 0*ones(size(transparent)), 1*ones(size(transparent)));

data.hcolor_saggital = imagesc(colour); 
hold off
set(data.hcolor_saggital, 'AlphaData', transparent*data.color_strength)
%\----------------------------------------------------------------%
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Saggital View: xslice = '},num2str(posx), {' out of '} ,num2str(data.sizes(1))))
% WARNING: I TRANSPOSED THE IMAGE, SO X AND Y ARE FLIPPED!!!!!

setMyData(data)

% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_show_slices_sag_overlay_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data = getMyData();
data.posx = get(handles.slider3,'Value');
posx=round(data.sizes(1)*data.posx);
if posx == 0; posx = 1; end
% Display Image Now
set(gcf,'CurrentAxes',handles.axes3)
if(~isempty(data.clims))
    imagesc(squeeze(data.voxelvolume(posx,:,:))',[data.clims(1) data.clims(2)])
else
    imagesc(squeeze(data.voxelvolume(posx,:,:))')
end
%-----------Include this part for the transparent overlay---------%
hold on
transparent = scale_uint8(  squeeze(data.voxelvolume2(posx,:,:))  )';  

% set transparanecy color in rgb.
colour = cat(3, 1*ones(size(transparent)), 0*ones(size(transparent)), 1*ones(size(transparent)));

data.hcolor_saggital = imagesc(colour); 
hold off
set(data.hcolor_saggital, 'AlphaData', transparent*data.color_strength)
%\----------------------------------------------------------------%
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
% colormap(hot)
title(strcat({'Saggital View: xslice = '},num2str(posx), {' out of '} ,num2str(data.sizes(1))))
setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function data=getMyData()
data=getappdata(gcf,'showdata');
% data = getappdata(gca,'showdata'); same result, just need to be
% consistent with what handle was passed in 'setMyData'.

function setMyData(data)
setappdata(data.HandleWindow,'showdata',data);
% setappdata(gca,'showdata',data); same result
