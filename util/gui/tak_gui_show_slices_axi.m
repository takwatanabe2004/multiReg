function varargout = tak_gui_show_slices_axi(varargin)
%% This function displays the slices of a 3d volume from the 3-coordinates...i.e. the 
%  axial (xy), coronal (xz), and saggital (yz) view.
%  varargin{1} = image volume to display slices
%  varargin{2} = clim = [lo hi]
%  varargin{3} = slider spacing setting (#s of index to slide) = [smallz bigz].
%  varargin{4} = [data.cbar data.cmap data.fill]

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tak_gui_show_slices_axi_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_show_slices_axi_OutputFcn, ...
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

% --- Executes just before tak_gui_show_slices_axi is made visible.
function tak_gui_show_slices_axi_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tak_gui_show_slices_axi
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% - Default value for setting
data.cbar = 0; % if cbara ==1, display colorbar on top
data.cmap = 0; % if cmap == 1, show image with color
data.fill = 1; % if fill == 1, fill axis (ie axis('normal'))

if(size(varargin,2)==4) && (~isempty(varargin{4}))
    data.cbar = varargin{4}.cbar; 
    data.cmap = varargin{4}.cmap; 
    data.fill = varargin{4}.fill; 
end

data.HandleWindow = gcf;
% set(data.HandleWindow,'Renderer','opengl');
data.voxelvolume = varargin{1};
data.sizes = size(data.voxelvolume);
data.posx = 0.25; data.posy = 0.5; data.posz = 0.5;

data.clims = [];    
if(size(varargin,2)>=2) && (~isempty(varargin{2}))
    data.clims = varargin{2};    
end

if(size(varargin,2)>=3) && (~isempty(varargin{3}))
% Set slider value range accordingly
% The order is in a hokie way since i wanted the input to be in the form  
% [smallx bigx smally bigy smallz bigz]
%     data.spacing3_small = varargin{3}(1)/data.sizes(1);
%     data.spacing3_big   = varargin{3}(2)/data.sizes(1);
%     data.spacing2_small = varargin{3}(3)/data.sizes(2);
%     data.spacing2_big   = varargin{3}(4)/data.sizes(2);
    data.spacing1_small = varargin{3}(1)/data.sizes(3);
    data.spacing1_big   = varargin{3}(2)/data.sizes(3);    
    set(handles.slider1,'SliderStep',[data.spacing1_small data.spacing1_big])
%     set(handles.slider2,'SliderStep',[data.spacing2_small data.spacing2_big])
%     set(handles.slider3,'SliderStep',[data.spacing3_small data.spacing3_big])
else
    data.spacing1_small = 1/data.sizes(3);
    data.spacing1_big   = 10/data.sizes(3);    
    set(handles.slider1,'SliderStep',[data.spacing1_small data.spacing1_big])    
end



% Now set the handles and figure according to the slider location
set(handles.slider1,'value',data.posx);
% set(handles.slider2,'value',data.posy);
% set(handles.slider3,'value',data.posz);

posx=round(data.sizes(1)*data.posx);
posy=round(data.sizes(2)*data.posy);
posz=round(data.sizes(3)*data.posz);

% Display axial view
set(gcf,'CurrentAxes',handles.axes1)
if(~isempty(data.clims))
    imagesc(data.voxelvolume(:,:,posz)',[data.clims(1) data.clims(2)])
else
    imagesc(data.voxelvolume(:,:,posz)')
end
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Axial View: zslice = '},num2str(posz), {' out of '} ,num2str(data.sizes(3))))

% % Display Coronal view
% set(gcf,'CurrentAxes',handles.axes2)
% if(~isempty(data.clims))
%     imagesc(squeeze(data.voxelvolume(:,posy,:))',[data.clims(1) data.clims(2)])
% else    
%     imagesc(squeeze(data.voxelvolume(:,posy,:))')
% end
% axis('xy','off','image')
% if data.cbar == 1; colorbar('location','northoutside'); end
% if data.cmap == 0; colormap('gray'); end   
% if data.fill == 1; axis('normal'); end   
% title(strcat({'Coronal View: yslice = '},num2str(posy), {' out of '} ,num2str(data.sizes(2))))
% 
% % Display Saggital view
% set(gcf,'CurrentAxes',handles.axes3)
% if(~isempty(data.clims))
%     imagesc(squeeze(data.voxelvolume(posx,:,:))',[data.clims(1) data.clims(2)])
% else
%     imagesc(squeeze(data.voxelvolume(posx,:,:))')
% end
% axis('xy','off','image')
% if data.cbar == 1; colorbar('location','northoutside'); end
% if data.cmap == 0; colormap('gray'); end   
% if data.fill == 1; axis('normal'); end   
% title(strcat({'Saggital View: xslice = '},num2str(posx), {' out of '} ,num2str(data.sizes(1))))
% WARNING: I TRANSPOSED THE IMAGE, SO X AND Y ARE FLIPPED!!!!!

setMyData(data)

% UIWAIT makes tak_gui_show_slices_axi wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_show_slices_axi_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
data = getMyData();
data.posz = get(handles.slider1,'Value');
posz=round(data.sizes(3)*data.posz);
if posz == 0; posz = 1; end
% Display Image Now
set(gcf,'CurrentAxes',handles.axes1)
if(~isempty(data.clims))
    imagesc(data.voxelvolume(:,:,posz)',[data.clims(1) data.clims(2)])
else
    imagesc(data.voxelvolume(:,:,posz)')
end
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end  
if data.fill == 1; axis('normal'); end   
title(strcat({'Axial View: zslice = '},num2str(posz), {' out of '} ,num2str(data.sizes(3))))
setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% % --- Executes on slider movement.
% function slider2_Callback(hObject, eventdata, handles)
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% data = getMyData();
% data.posy = get(handles.slider2,'Value');
% posy=round(data.sizes(2)*data.posy);
% if posy == 0; posy = 1; end
% 
% % Display Image Now
% set(gcf,'CurrentAxes',handles.axes2)
% if(~isempty(data.clims))
%     imagesc(squeeze(data.voxelvolume(:,posy,:))',[data.clims(1) data.clims(2)])
% else    
%     imagesc(squeeze(data.voxelvolume(:,posy,:))')
% end
% axis('xy','off','image')
% if data.cbar == 1; colorbar('location','northoutside'); end
% if data.cmap == 0; colormap('gray'); end   
% if data.fill == 1; axis('normal'); end   
% title(strcat({'Coronal View: yslice = '},num2str(posy), {' out of '} ,num2str(data.sizes(2))))
% setMyData(data)
% 
% % --- Executes during object creation, after setting all properties.
% function slider2_CreateFcn(hObject, eventdata, handles)
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end
% 
% 
% % --- Executes on slider movement.
% function slider3_Callback(hObject, eventdata, handles)
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% data = getMyData();
% data.posx = get(handles.slider3,'Value');
% posx=round(data.sizes(1)*data.posx);
% if posx == 0; posx = 1; end
% % Display Image Now
% set(gcf,'CurrentAxes',handles.axes3)
% if(~isempty(data.clims))
%     imagesc(squeeze(data.voxelvolume(posx,:,:))',[data.clims(1) data.clims(2)])
% else
%     imagesc(squeeze(data.voxelvolume(posx,:,:))')
% end
% axis('xy','off','image')
% if data.cbar == 1; colorbar('location','northoutside'); end
% if data.cmap == 0; colormap('gray'); end   
% if data.fill == 1; axis('normal'); end   
% % colormap(hot)
% title(strcat({'Saggital View: xslice = '},num2str(posx), {' out of '} ,num2str(data.sizes(1))))
% setMyData(data)
% 
% % --- Executes during object creation, after setting all properties.
% function slider3_CreateFcn(hObject, eventdata, handles)
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end

function data=getMyData()
data=getappdata(gcf,'showdata');
% data = getappdata(gca,'showdata'); same result, just need to be
% consistent with what handle was passed in 'setMyData'.

function setMyData(data)
setappdata(data.HandleWindow,'showdata',data);
% setappdata(gca,'showdata',data); same result
