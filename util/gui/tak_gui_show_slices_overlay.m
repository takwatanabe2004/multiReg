function varargout = tak_gui_show_slices_overlay(varargin)
%% This function displays the slices of a 3d volume from the 3-coordinates...i.e. the 
%  axial (xy), coronal (xz), and saggital (yz) view.
%  varargin{1} = image volume to display slices
%  varargin{2} = image volume to overlay
%  varargin{3} = overlay color transparency strength (default 0.25)
%  varargin{4} = clim = [lo hi] (defautl [])
%  varargin{5} = initial slider position [posz, posy, posx]
%                => (this order is to correspond to the slider dimension direction)
%  varargin{6} = slider spacing setting (#s of index to slide) = [smallx bigx smally bigy smallz bigz].
%  varargin{7} = structure data (data.cbar, data.cmap, data.fill data.sliderpos)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tak_gui_show_slices_overlay_OpeningFcn, ...
                   'gui_OutputFcn',  @tak_gui_show_slices_overlay_OutputFcn, ...
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

% --- Executes just before tak_gui_show_slices2 is made visible.
function tak_gui_show_slices_overlay_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tak_gui_show_slices2
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%---------- Process input data first ----------------------%
% main volume to display
data.voxelvolume = varargin{1};
data.sizes = size(data.voxelvolume);

% volume to overlay
data.voxelvolume2 = varargin{2};

if (size(varargin,2)>=3) && (~isempty(varargin{3}))
    data.color_strength = varargin{3};
else
    data.color_strength = 0.25;
end

if(size(varargin,2)>=4) && (~isempty(varargin{4}))
    data.clims = varargin{4};    
else
    data.clims = [];    
end

if(size(varargin,2)>=5) && (~isempty(varargin{5}))
    data.posz = varargin{5}(1);
    data.posy = varargin{5}(2);
    data.posx = varargin{5}(3);
else
    data.posx = 0.5; data.posy = 0.5; data.posz = 0.5;
end

if(size(varargin,2)>=6) && (~isempty(varargin{6}))
% Set slider value range accordingly
% The order is in a hokie way since i wanted the input to be in the form  
% [smallx bigx smally bigy smallz bigz]
    data.spacing3_small = varargin{6}(1)/data.sizes(1);
    data.spacing3_big   = varargin{6}(2)/data.sizes(1);
    data.spacing2_small = varargin{6}(3)/data.sizes(2);
    data.spacing2_big   = varargin{6}(4)/data.sizes(2);
    data.spacing1_small = varargin{6}(5)/data.sizes(3);
    data.spacing1_big   = varargin{6}(6)/data.sizes(3);    
    set(handles.slider1,'SliderStep',[data.spacing1_small data.spacing1_big])
    set(handles.slider2,'SliderStep',[data.spacing2_small data.spacing2_big])
    set(handles.slider3,'SliderStep',[data.spacing3_small data.spacing3_big])
else % default
    data.spacing3_small = 1/data.sizes(1);
    data.spacing3_big   = 10/data.sizes(1);
    data.spacing2_small = 1/data.sizes(2);
    data.spacing2_big   = 10/data.sizes(2);
    data.spacing1_small = 1/data.sizes(3);
    data.spacing1_big   = 10/data.sizes(3);    
    set(handles.slider1,'SliderStep',[data.spacing1_small data.spacing1_big])
    set(handles.slider2,'SliderStep',[data.spacing2_small data.spacing2_big])
    set(handles.slider3,'SliderStep',[data.spacing3_small data.spacing3_big])
end

if(size(varargin,2)>=7) && (~isempty(varargin{7}))
    data.cbar = varargin{7}.cbar; 
    data.cmap = varargin{7}.cmap; 
    data.fill = varargin{7}.fill; 
else
    data.cbar = 0; % if cbara == 1, display colorbar on top
    data.cmap = 0; % if cmap == 1, show image with color
    data.fill = 1; % if fill == 1, fill axis (ie axis('normal'))
end
%\-------------------------------------------------%
data.HandleWindow = gcf;
% set(data.HandleWindow,'Renderer','opengl');

data.xline = linspace(0,data.sizes(1),2);
data.xones = ones(size(data.xline));
data.yline = linspace(0,data.sizes(2),2);
data.yones = ones(size(data.yline));
data.zline = linspace(0,data.sizes(3),2);
data.zones = ones(size(data.zline));

% Now set the handles and figure according to the slider location
set(handles.slider1,'value',data.posz);
set(handles.slider2,'value',data.posy);
set(handles.slider3,'value',data.posx);

posx=round(data.sizes(1)*data.posx);
posy=round(data.sizes(2)*data.posy);
posz=round(data.sizes(3)*data.posz);

% Display axial view
set(gcf,'CurrentAxes',handles.axes1)
if(~isempty(data.clims))
    data.h_axi = imagesc(data.voxelvolume(:,:,posz)',[data.clims(1) data.clims(2)]);
else
    data.h_axi = imagesc(data.voxelvolume(:,:,posz)');
end
%-----------Include this part for the transparent overlay---------%
hold on
transparent = scale_uint8(data.voxelvolume2(:,:,posz))';  

% set transparanecy color in rgb.
colour = cat(3, 1*ones(size(transparent)), 0*ones(size(transparent)), 1*ones(size(transparent)));

data.hcolor_axial = imagesc(colour); 
hold off
set(data.hcolor_axial, 'AlphaData', transparent*data.color_strength)
%\----------------------------------------------------------------%
hold on
data.h_xl_yp = line(data.xline,posy*data.xones,'Color','b');
data.h_xp_yl = line(posx*data.yones,data.yline,'Color','r');
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Axial View: zslice = '},num2str(posz), {' out of '} ,num2str(data.sizes(3))))

% Display Coronal view
set(gcf,'CurrentAxes',handles.axes2)
if(~isempty(data.clims))
    data.h_cor = imagesc(squeeze(data.voxelvolume(:,posy,:))',[data.clims(1) data.clims(2)]);
else    
    data.h_cor = imagesc(squeeze(data.voxelvolume(:,posy,:))');
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
hold on
data.h_xl_zp = line(data.xline,posz*data.xones,'Color','b');
data.h_xp_zl = line(posx*data.zones,data.zline,'Color','r');
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Coronal View: yslice = '},num2str(posy), {' out of '} ,num2str(data.sizes(2))))

% Display Saggital view
set(gcf,'CurrentAxes',handles.axes3)
if(~isempty(data.clims))
    data.h_sag = imagesc(squeeze(data.voxelvolume(posx,:,:))',[data.clims(1) data.clims(2)]);
else
    data.h_sag = imagesc(squeeze(data.voxelvolume(posx,:,:))');
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
hold on
data.h_yl_zp = line(data.yline,posz*data.yones,'Color','b');
data.h_yp_zl = line(posy*data.zones,data.zline,'Color','r');
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end   
if data.fill == 1; axis('normal'); end   
title(strcat({'Saggital View: xslice = '},num2str(posx), {' out of '} ,num2str(data.sizes(1))))
% WARNING: I TRANSPOSED THE IMAGE, SO X AND Y ARE FLIPPED!!!!!

setMyData(data)

% UIWAIT makes tak_gui_show_slices2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = tak_gui_show_slices_overlay_OutputFcn(hObject, eventdata, handles) 
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
    data.h_axi = imagesc(data.voxelvolume(:,:,posz)',[data.clims(1) data.clims(2)]);
else
    data.h_axi = imagesc(data.voxelvolume(:,:,posz)');
end
%-----------Include this part for the transparent overlay---------%
hold on
transparent = scale_uint8(data.voxelvolume2(:,:,posz))';  

% set transparanecy color in rgb.
colour = cat(3, 1*ones(size(transparent)), 0*ones(size(transparent)), 1*ones(size(transparent)));

data.hcolor_axial = imagesc(colour); 
hold off
set(data.hcolor_axial, 'AlphaData', transparent*data.color_strength)
%\----------Include this part for the transparent overlay---------%
axis('xy','off','image')
if data.cbar == 1; colorbar('location','northoutside'); end
if data.cmap == 0; colormap('gray'); end  
if data.fill == 1; axis('normal'); end   
title(strcat({'Axial View: zslice = '},num2str(posz), {' out of '} ,num2str(data.sizes(3))))

%---------------Modify other lines-------------------%
data.posy = get(handles.slider2,'Value');
posy=round(data.sizes(2)*data.posy);
if posy == 0; posy = 1; end

data.posx = get(handles.slider3,'Value');
posx=round(data.sizes(1)*data.posx);
if posx == 0; posx = 1; end

hold on
data.h_xl_yp = line(data.xline,posy*data.xones,'Color','b');
data.h_xp_yl = line(posx*data.yones,data.yline,'Color','r');

% fix coronal image
set(gcf,'CurrentAxes',handles.axes2)
cla(data.h_xl_zp,[data.h_cor; data.h_xp_zl;data.hcolor_coronal])
hold on
data.h_xl_zp = line(data.xline,posz*data.xones,'Color','b');

% fix sagital image
set(gcf,'CurrentAxes',handles.axes3)
cla(data.h_yl_zp,[data.h_sag;data.h_yp_zl;data.hcolor_saggital])
hold on
data.h_yl_zp = line(data.yline,posz*data.yones,'Color','b');

setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


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
    data.h_cor = imagesc(squeeze(data.voxelvolume(:,posy,:))',[data.clims(1) data.clims(2)]);
else    
    data.h_cor = imagesc(squeeze(data.voxelvolume(:,posy,:))');
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

%---------------Modify other lines-------------------%
data.posz = get(handles.slider1,'Value');
posz=round(data.sizes(3)*data.posz);
if posz == 0; posz = 1; end

data.posx = get(handles.slider3,'Value');
posx=round(data.sizes(1)*data.posx);
if posx == 0; posx = 1; end

hold on
data.h_xl_zp = line(data.xline,posz*data.xones,'Color','b');
data.h_xp_zl = line(posx*data.zones,data.zline,'Color','r');

% fix axial image
set(gcf,'CurrentAxes',handles.axes1)
cla(data.h_xl_yp,[data.h_axi; data.h_xp_yl;data.hcolor_axial])
hold on
data.h_xl_yp = line(data.xline,posy*data.xones,'Color','b');

% fix sagital image
set(gcf,'CurrentAxes',handles.axes3)
cla(data.h_yp_zl,[data.h_sag;data.h_yl_zp;data.hcolor_saggital])
hold on
data.h_yp_zl = line(posy*data.zones,data.zline,'Color','r');

setMyData(data)

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


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
    data.h_sag = imagesc(squeeze(data.voxelvolume(posx,:,:))',[data.clims(1) data.clims(2)]);
else
    data.h_sag = imagesc(squeeze(data.voxelvolume(posx,:,:))');
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

%---------------Modify other lines-------------------%
data.posz = get(handles.slider1,'Value');
posz=round(data.sizes(3)*data.posz);
if posz == 0; posz = 1; end

data.posy = get(handles.slider2,'Value');
posy=round(data.sizes(2)*data.posy);
if posy == 0; posy = 1; end

hold on
data.h_yl_zp = line(data.yline,posz*data.yones,'Color','b');
data.h_yp_zl = line(posy*data.zones,data.zline,'Color','r');

% fix axial image
set(gcf,'CurrentAxes',handles.axes1)
cla(data.h_xp_yl,[data.h_axi; data.h_xl_yp;data.hcolor_axial])
hold on
data.h_xp_yl = line(posx*data.yones,data.yline,'Color','r');

% fix coronal image
set(gcf,'CurrentAxes',handles.axes2)
cla(data.h_xp_zl,[data.h_cor;data.h_xl_zp;data.hcolor_coronal])
hold on
data.h_xp_zl = line(posx*data.zones,data.zline,'Color','r');

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
