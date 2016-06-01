function varargout = CompMapGenerator(varargin)
% COMPMAPGENERATOR MATLAB code for CompMapGenerator.fig
%      COMPMAPGENERATOR, by itself, creates a new COMPMAPGENERATOR or raises the existing
%      singleton*.
%
%      H = COMPMAPGENERATOR returns the handle to a new COMPMAPGENERATOR or the handle to
%      the existing singleton*.
%
%      COMPMAPGENERATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPMAPGENERATOR.M with the given input arguments.
%
%      COMPMAPGENERATOR('Property','Value',...) creates a new COMPMAPGENERATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CompMapGenerator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CompMapGenerator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CompMapGenerator

% Last Modified by GUIDE v2.5 13-Apr-2015 23:04:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CompMapGenerator_OpeningFcn, ...
                   'gui_OutputFcn',  @CompMapGenerator_OutputFcn, ...
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


% --- Executes just before CompMapGenerator is made visible.
function CompMapGenerator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CompMapGenerator (see VARARGIN)

% Choose default command line output for CompMapGenerator
handles.output = hObject;
handles.SAE_map_compOrg = varargin{1};
handles.SAE_map_compNew = varargin{1};
if length(varargin) > 1
    handles.compRef = varargin{2};
else
    handles.compRef = [];
end;
    
if isempty(handles.compRef)
    DrawComp(handles.SAE_map_compOrg,handles.CompMap,handles.compRef(:,1),handles.compRef(:,2),1e5,298,[]);
else
    DrawComp(handles.SAE_map_compOrg,handles.CompMap,handles.compRef(:,1),handles.compRef(:,2),1e5,298,[]);
end;
RPMRep = num2cell(unique(handles.SAE_map_compOrg.RPM));
set(handles.RPMLB,'String',RPMRep);
step = str2double(get(handles.RPMStepEdit,'String'));
RPMNew = [RPMRep{1}; ...
    (ceil(RPMRep{1}/step)*step:step:RPMRep{end})';RPMRep{end}];
set(handles.RPMNewLB,'String',RPMNew);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CompMapGenerator wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CompMapGenerator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.SAE_map_compNew;
close(handles.figure1);

% --- Executes on selection change in RPMLB.
function RPMLB_Callback(hObject, eventdata, handles)
% hObject    handle to RPMLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RPMLB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RPMLB


% --- Executes during object creation, after setting all properties.
function RPMLB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RPMLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in regenSAEPB.
function regenSAEPB_Callback(hObject, eventdata, handles)
% hObject    handle to regenSAEPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RPMNew = cellfun(@str2double,cellstr(get(handles.RPMNewLB,'String')));
[handles.SAE_map_compNew] ...
    = produceNewSAEComp(RPMNew,handles.SAE_map_compOrg);
if isempty(handles.compRef)
    DrawComp(handles.SAE_map_compNew,handles.CompMap,handles.compRef(:,1),handles.compRef(:,2),1e5,298,[]);
else
    DrawComp(handles.SAE_map_compNew,handles.CompMap,handles.compRef(:,1),handles.compRef(:,2),1e5,298,[]);
end;
handles.SAE_map_compOrg = handles.SAE_map_compNew;
guidata(hObject,handles);

% --- Executes on button press in exportMapPB.
function exportMapPB_Callback(hObject, eventdata, handles)
% hObject    handle to exportMapPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);



function RPMStepEdit_Callback(hObject, eventdata, handles)
% hObject    handle to RPMStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
step = str2double(get(hObject,'String'));
RPMRep = cellfun(@str2double,cellstr(get(handles.RPMNewLB,'String')));
RPMNew = [RPMRep(1); ...
    (ceil(RPMRep(1)/step)*step:step:RPMRep(end))';RPMRep(end)];
set(handles.RPMNewLB,'String',RPMNew);

% Hints: get(hObject,'String') returns contents of RPMStepEdit as text
%        str2double(get(hObject,'String')) returns contents of RPMStepEdit as a double


% --- Executes during object creation, after setting all properties.
function RPMStepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RPMStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RPMNewLB.
function RPMNewLB_Callback(hObject, eventdata, handles)
% hObject    handle to RPMNewLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RPMNewLB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RPMNewLB


% --- Executes during object creation, after setting all properties.
function RPMNewLB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RPMNewLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mDotResizeRatioEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mDotResizeRatioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mDotResizeRatioEdit as text
%        str2double(get(hObject,'String')) returns contents of mDotResizeRatioEdit as a double


% --- Executes during object creation, after setting all properties.
function mDotResizeRatioEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mDotResizeRatioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PRResizeRatioEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PRResizeRatioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PRResizeRatioEdit as text
%        str2double(get(hObject,'String')) returns contents of PRResizeRatioEdit as a double


% --- Executes during object creation, after setting all properties.
function PRResizeRatioEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PRResizeRatioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ResizePB.
function ResizePB_Callback(hObject, eventdata, handles)
% hObject    handle to ResizePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mDotResizeRatio = str2num(get(handles.mDotResizeRatioEdit,'String'));
PRResizeRatio = str2num(get(handles.PRResizeRatioEdit,'String'));
handles.SAE_map_compNew.PR = handles.SAE_map_compOrg.PR ...
    *PRResizeRatio/100;
handles.SAE_map_compNew.PR_chock = handles.SAE_map_compOrg.PR_chock ...
    *PRResizeRatio/100;
handles.SAE_map_compNew.PR_surge = handles.SAE_map_compOrg.PR_surge ...
    *PRResizeRatio/100;

handles.SAE_map_compNew.m_dot = handles.SAE_map_compOrg.m_dot ...
    *mDotResizeRatio/100;
handles.SAE_map_compNew.m_dot_surge = handles.SAE_map_compOrg.m_dot_surge ...
    *mDotResizeRatio/100;
handles.SAE_map_compNew.m_dot_chock = handles.SAE_map_compOrg.m_dot_chock ...
    *mDotResizeRatio/100;

RPMNew = cellfun(@str2double,cellstr(get(handles.RPMNewLB,'String')));
[handles.SAE_map_compNew] ...
    = produceNewSAEComp(RPMNew,handles.SAE_map_compNew);
if isempty(handles.compRef)
    DrawComp(handles.SAE_map_compNew,handles.CompMap,handles.compRef(:,1),handles.compRef(:,2),1e5,298,[]);
else
    DrawComp(handles.SAE_map_compNew,handles.CompMap,handles.compRef(:,1),handles.compRef(:,2),1e5,298,[]);
end;

guidata(hObject,handles);
