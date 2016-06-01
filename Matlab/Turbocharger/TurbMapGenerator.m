function varargout = TurbMapGenerator(varargin)
% TURBMAPGENERATOR MATLAB code for TurbMapGenerator.fig
%      TURBMAPGENERATOR, by itself, creates a new TURBMAPGENERATOR or raises the existing
%      singleton*.
%
%      H = TURBMAPGENERATOR returns the handle to a new TURBMAPGENERATOR or the handle to
%      the existing singleton*.
%
%      TURBMAPGENERATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TURBMAPGENERATOR.M with the given input arguments.
%
%      TURBMAPGENERATOR('Property','Value',...) creates a new TURBMAPGENERATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TurbMapGenerator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TurbMapGenerator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TurbMapGenerator

% Last Modified by GUIDE v2.5 09-Jun-2015 11:46:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TurbMapGenerator_OpeningFcn, ...
                   'gui_OutputFcn',  @TurbMapGenerator_OutputFcn, ...
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


% --- Executes just before TurbMapGenerator is made visible.
function TurbMapGenerator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TurbMapGenerator (see VARARGIN)

% Choose default command line output for TurbMapGenerator
handles.output = hObject;
handles.SAE_map_turbOrg = varargin{1};
handles.SAE_map_turbNew = varargin{1};
handles.turbRef = varargin{2};
if isempty(handles.turbRef)
    DrawTurb(handles.SAE_map_turbOrg,[],[],[],handles.CompMap);
else
    DrawTurb(handles.SAE_map_turbOrg,handles.turbRef(:,1),...
        handles.turbRef(:,2),handles.turbRef(:,3),handles.CompMap);
end;
RPMRep = num2cell(unique(handles.SAE_map_turbOrg.RPM));
set(handles.RPMLB,'String',RPMRep);
step = str2double(get(handles.RPMStepEdit,'String'));
RPMNew = [RPMRep{1}; ...
    (ceil(RPMRep{1}/step)*step:step:RPMRep{end})';RPMRep{end}];
set(handles.RPMNewLB,'String',RPMNew);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TurbMapGenerator wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TurbMapGenerator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.SAE_map_turbNew;
%close(handles.figure1);

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
[handles.SAE_map_turbNew] ...
    = produceNewSAEturb(RPMNew,handles.SAE_map_turbOrg);
if isempty(handles.turbRef)
    DrawTurb(handles.SAE_map_turbNew,[],[],[],handles.CompMap);
else
    DrawTurb(handles.SAE_map_turbNew,handles.turbRef(:,1),...
        handles.turbRef(:,2),handles.turbRef(:,3),handles.CompMap);
end;
handles.SAE_map_turbOrg = handles.SAE_map_turbNew;
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
mDotResizeShift = str2num(get(handles.mDotResizeShiftEdit,'String'));
mDotResizeShift = (mDotResizeShift/100-1)*handles.SAE_map_turbOrg.m_dot(1);
PRResizeRatio = str2num(get(handles.PRResizeRatioEdit,'String'));
PRResizeShift = str2num(get(handles.PRResizeShiftEdit,'String'));
PRResizeShift = (PRResizeShift/100-1)*handles.SAE_map_turbOrg.PR(1);

handles.SAE_map_turbNew.PR = handles.SAE_map_turbOrg.PR ...
    *PRResizeRatio/100 + PRResizeShift;

handles.SAE_map_turbNew.m_dot = handles.SAE_map_turbOrg.m_dot ...
    *mDotResizeRatio/100 + mDotResizeShift;

%RPMNew = cellfun(@str2double,cellstr(get(handles.RPMNewLB,'String')));
%[handles.SAE_map_turbNew] ...
%    = produceNewSAEComp(RPMNew,handles.SAE_map_turbNew);
if isempty(handles.turbRef)
    DrawTurb(handles.SAE_map_turbNew,[],[],[],handles.CompMap);
else
    DrawTurb(handles.SAE_map_turbNew,handles.turbRef(:,1),...
        handles.turbRef(:,2),handles.turbRef(:,3),handles.CompMap);
end;

guidata(hObject,handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function mDotResizeShiftEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mDotResizeShiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mDotResizeShiftEdit as text
%        str2double(get(hObject,'String')) returns contents of mDotResizeShiftEdit as a double


% --- Executes during object creation, after setting all properties.
function mDotResizeShiftEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mDotResizeShiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PRResizeShiftEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PRResizeShiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PRResizeShiftEdit as text
%        str2double(get(hObject,'String')) returns contents of PRResizeShiftEdit as a double


% --- Executes during object creation, after setting all properties.
function PRResizeShiftEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PRResizeShiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
