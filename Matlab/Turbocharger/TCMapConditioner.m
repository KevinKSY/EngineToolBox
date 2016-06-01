function varargout = TCMapConditioner(varargin)
% TCMAPCONDITIONER MATLAB code for TCMapConditioner.fig
%      TCMAPCONDITIONER, by itself, creates a new TCMAPCONDITIONER or raises the existing
%      singleton*.
%
%      H = TCMAPCONDITIONER returns the handle to a new TCMAPCONDITIONER or the handle to
%      the existing singleton*.
%
%      TCMAPCONDITIONER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TCMAPCONDITIONER.M with the given input arguments.
%
%      TCMAPCONDITIONER('Property','Value',...) creates a new TCMAPCONDITIONER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TCMapConditioner_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TCMapConditioner_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TCMapConditioner

% Last Modified by GUIDE v2.5 10-Apr-2015 15:57:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TCMapConditioner_OpeningFcn, ...
                   'gui_OutputFcn',  @TCMapConditioner_OutputFcn, ...
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


% --- Executes just before TCMapConditioner is made visible.
function TCMapConditioner_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TCMapConditioner (see VARARGIN)

% Choose default command line output for TCMapConditioner
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TCMapConditioner wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TCMapConditioner_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function CompMapFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CompMapFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CompMapFileEdit as text
%        str2double(get(hObject,'String')) returns contents of CompMapFileEdit as a double


% --- Executes during object creation, after setting all properties.
function CompMapFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CompMapFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindCompFilePB.
function FindCompFilePB_Callback(hObject, eventdata, handles)
% hObject    handle to FindCompFilePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName, pathName] = uigetfile({'*.mat'});
set(handles.CompMapFileEdit,'String',fileName);
handles.compFile = {fileName,pathName};
guidata(hObject, handles);


function TurbMapFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TurbMapFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TurbMapFileEdit as text
%        str2double(get(hObject,'String')) returns contents of TurbMapFileEdit as a double


% --- Executes during object creation, after setting all properties.
function TurbMapFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TurbMapFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindTurbFilePB.
function FindTurbFilePB_Callback(hObject, eventdata, handles)
% hObject    handle to FindTurbFilePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName, pathName] = uigetfile({'*.mat'});
set(handles.TurbMapFileEdit,'String',fileName);
handles.turbFile = {fileName,pathName};
guidata(hObject, handles);

% --- Executes on button press in MapLoadPB.
function MapLoadPB_Callback(hObject, eventdata, handles)
% hObject    handle to MapLoadPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load([handles.compFile{2} handles.compFile{1}]);
load([handles.turbFile{2} handles.turbFile{1}]);
handles.SAE_map_comp = SAE_map_comp;
handles.SAE_map_turb = SAE_map_turb;
cla(handles.CompMap);
cla(handles.TurbMap);
line(SAE_map_comp.m_dot,SAE_map_comp.PR,SAE_map_comp.eff,...
    'marker','.','markersize',4,'linestyle','none','parent',...
    handles.CompMap);
line(SAE_map_turb.PR,SAE_map_turb.m_dot,SAE_map_turb.eff,...
    'marker','.','markersize',4,'linestyle','none','parent',...
    handles.TurbMap);
handles.RPMComp_rep = num2cell(unique(SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;
set(handles.RPMComplistbox,'String',handles.RPMComp_rep);
set(handles.noPtCompEdit,'String',length(handles.compRPMIdx{1}));
set(handles.RPMTurblistbox,'String',handles.RPMTurb_rep);
set(handles.noPtTurbEdit,'String',length(handles.turbRPMIdx{1}));
cla(handles.CompCondInterp);
cla(handles.TurbCondInterp);
plot(handles.CompCondInterp,SAE_map_comp.PR(handles.compRPMIdx{1}),...
    handles.SAE_map_comp.m_dot(handles.compRPMIdx{1}),'*r');
plot(handles.TurbCondInterp,SAE_map_turb.PR(handles.turbRPMIdx{1}),...
    handles.SAE_map_turb.m_dot(handles.turbRPMIdx{1}),'*r');
guidata(hObject,handles);


% --- Executes on selection change in RPMTurblistbox.
function RPMTurblistbox_Callback(hObject, eventdata, handles)
% hObject    handle to RPMTurblistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;
IdxSelected = get(hObject,'Value');
RPMSelected = handles.RPMTurb_rep{IdxSelected};
set(handles.noPtTurbEdit,'String',length(handles.turbRPMIdx{IdxSelected}));
SB = get(get(handles.turbDataBG,'SelectedObject'),'tag');
cla(handles.TurbCondInterp);
switch SB
    case 'turbMassFlowRB'
        plot(handles.TurbCondInterp,...
            handles.SAE_map_turb.PR(handles.turbRPMIdx{IdxSelected}),...
            handles.SAE_map_turb.m_dot(handles.turbRPMIdx{IdxSelected})...
            ,'*r');    
    case 'turbEffRB'
        plot(handles.TurbCondInterp,...
            handles.SAE_map_turb.PR(handles.turbRPMIdx{IdxSelected}),...
            handles.SAE_map_turb.eff(handles.turbRPMIdx{IdxSelected})...
            ,'*r');    
end;






% Hints: contents = cellstr(get(hObject,'String')) returns RPMTurblistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RPMTurblistbox


% --- Executes during object creation, after setting all properties.
function RPMTurblistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RPMTurblistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noPtCompEdit_Callback(hObject, eventdata, handles)
% hObject    handle to noPtCompEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noPtCompEdit as text
%        str2double(get(hObject,'String')) returns contents of noPtCompEdit as a double


% --- Executes during object creation, after setting all properties.
function noPtCompEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noPtCompEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in turbInterpPB.
function turbInterpPB_Callback(hObject, eventdata, handles)
% hObject    handle to turbInterpPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;

IdxSelected = get(handles.RPMTurblistbox,'Value');
RPMSelected = handles.RPMTurb_rep{IdxSelected};
noPt = str2num(get(handles.noPtTurbEdit,'String'));
SB1 = get(get(handles.turbDataBG,'SelectedObject'),'tag');
SB2 = get(get(handles.turbInterpMethodBG,'SelectedObject'),'tag');
cla(handles.TurbCondInterp);
switch SB1
    case 'turbMassFlowRB'
        x = handles.SAE_map_turb.PR(handles.turbRPMIdx{IdxSelected});
        y = handles.SAE_map_turb.m_dot(handles.turbRPMIdx{IdxSelected});
        plot(handles.TurbCondInterp,x,y,'O'); 
        hold on
        switch SB2
            case 'turbLinearRB'
                interpMet = 'linear';
            case 'turbCubicRB'
                interpMet = 'cubic';
            case 'turbPCHPRB'
                interpMet = 'pchip';
        end;
        xNode = linspace(min(x),max(x),noPt);
        yInt = interp1(x,y,xNode',interpMet);
        line(xNode,yInt,'marker','.','markersize',4,'linestyle','none','parent',handles.TurbCondInterp);  
        handles.tempTurbRPM = RPMSelected*ones(noPt,1);
        handles.tempTurbMDot = yInt;
        handles.tempTurbPR = xNode';
    case 'turbEffRB'
        x = handles.SAE_map_turb.PR(handles.turbRPMIdx{IdxSelected});
        y = handles.SAE_map_turb.eff(handles.turbRPMIdx{IdxSelected});
        plot(handles.TurbCondInterp,x,y,'O'); 
        hold on
        switch SB2
            case 'turbLinearRB'
                interpMet = 'linear';
            case 'turbCubicRB'
                interpMet = 'cubic';
            case 'turbPCHPRB'
                interpMet = 'pchip';
        end;
        xNode = linspace(min(x),max(x),noPt);
        yInt = interp1(x,y,xNode',interpMet);
        line(xNode,yInt,'marker','.','markersize',4,'linestyle','none','parent',handles.TurbCondInterp);  
        handles.tempTurbRPM = RPMSelected*ones(noPt,1);
        handles.tempTurbEff = yInt;
        handles.tempTurbPR = xNode';
end;
guidata(hObject,handles);

% --- Executes on selection change in RPMComplistbox.
function RPMComplistbox_Callback(hObject, eventdata, handles)
% hObject    handle to RPMComplistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RPMComplistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RPMComplistbox
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;

IdxSelected = get(hObject,'Value');
RPMSelected = handles.RPMComp_rep{IdxSelected};
set(handles.noPtCompEdit,'String',length(handles.compRPMIdx{IdxSelected}));
SB = get(get(handles.compDataBG,'SelectedObject'),'tag');
cla(handles.CompCondInterp);
switch SB
    case 'compMassFlowRB'
        plot(handles.CompCondInterp,...
            handles.SAE_map_comp.PR(handles.compRPMIdx{IdxSelected}),...
            handles.SAE_map_comp.m_dot(handles.compRPMIdx{IdxSelected})...
            ,'*r');    
    case 'compEfficiencyRB'
        plot(handles.CompCondInterp,...
            handles.SAE_map_comp.PR(handles.compRPMIdx{IdxSelected}),...
            handles.SAE_map_comp.eff(handles.compRPMIdx{IdxSelected})...
            ,'*r');    
end;

% --- Executes during object creation, after setting all properties.
function RPMComplistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RPMComplistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TurbCondSavePB1.
function CompCondSavePB_Callback(hObject, eventdata, handles)
% hObject    handle to TurbCondSavePB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RPMTemp = handles.SAE_map_comp.RPM;
mDotTemp = handles.SAE_map_comp.m_dot;
PRTemp = handles.SAE_map_comp.PR;
EffTemp = handles.SAE_map_comp.eff;

idx = (RPMTemp == handles.tempRPM(1));
idx0 = find(idx(1:end-1) ~= idx(2:end));
idx1 = idx0(1)+1;
if (length(idx0) == 1)
    if (idx(1) == 1)
        mDotTemp = [handles.tempMDot;mDotTemp(idx1:end)];
        PRTemp = [handles.tempPR;PRTemp(idx1:end)];
        EffTemp = [handles.tempEff;EffTemp(idx1:end)];
        RPMTemp = [handles.tempRPM;RPMTemp(idx1:end)];
    else
        mDotTemp = [mDotTemp(1:idx1-1);handles.tempMDot];
        PRTemp = [PRTemp(1:idx1-1);handles.tempPR];
        EffTemp = [EffTemp(1:idx1-1);handles.tempEff];
        RPMTemp = [RPMTemp(1:idx1-1);handles.tempRPM];
    end;
else
    idx2 = idx0(2);
    mDotTemp = [mDotTemp(1:idx1-1);handles.tempMDot;mDotTemp(idx2+1:end)];
    PRTemp = [PRTemp(1:idx1-1);handles.tempPR;PRTemp(idx2+1:end)];
    EffTemp = [EffTemp(1:idx1-1);handles.tempEff;EffTemp(idx2+1:end)];
    RPMTemp = [RPMTemp(1:idx1-1);handles.tempRPM;RPMTemp(idx2+1:end)];
end;

handles.SAE_map_comp.m_dot = mDotTemp;
handles.SAE_map_comp.PR = PRTemp;
handles.SAE_map_comp.eff = EffTemp;
handles.SAE_map_comp.RPM = RPMTemp;

guidata(hObject,handles);



% --- Executes on button press in compInterpPB.
function compInterpPB_Callback(hObject, eventdata, handles)
% hObject    handle to compInterpPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;

IdxSelected = get(handles.RPMComplistbox,'Value');
RPMSelected = handles.RPMComp_rep{IdxSelected};
noPt = str2num(get(handles.noPtCompEdit,'String'));
SB1 = get(get(handles.compDataBG,'SelectedObject'),'tag');
SB2 = get(get(handles.compInterpMethodBG,'SelectedObject'),'tag');
cla(handles.CompCondInterp);
switch SB1
    case 'compMassFlowRB'
        x = handles.SAE_map_comp.PR(handles.compRPMIdx{IdxSelected});
        y = handles.SAE_map_comp.m_dot(handles.compRPMIdx{IdxSelected});
        plot(handles.CompCondInterp,x,y,'O'); 
        hold on
        switch SB2
            case 'compLinearRB'
                interpMet = 'linear';
            case 'compCubicRB'
                interpMet = 'cubic';
            case 'compPCHPRB'
                interpMet = 'pchip';
        end;
        xNode = linspace(min(x),max(x),noPt);
        yInt = interp1(x,y,xNode',interpMet);
        line(xNode,yInt,'marker','.','markersize',4,'linestyle','none','parent',handles.CompCondInterp);  
        handles.tempRPM = RPMSelected*ones(noPt,1);
        handles.tempMDot = yInt;
        handles.tempPR = xNode';
    case 'compEfficiencyRB'
        x = handles.SAE_map_comp.PR(handles.compRPMIdx{IdxSelected});
        y = handles.SAE_map_comp.eff(handles.compRPMIdx{IdxSelected});
        plot(handles.CompCondInterp,x,y,'O'); 
        hold on
        switch SB2
            case 'compLinearRB'
                interpMet = 'linear';
            case 'compCubicRB'
                interpMet = 'cubic';
            case 'compPCHPRB'
                interpMet = 'pchip';
        end;
        xNode = linspace(min(x),max(x),noPt);
        yInt = interp1(x,y,xNode',interpMet);
        line(xNode,yInt,'marker','.','markersize',4,'linestyle','none','parent',handles.CompCondInterp);  
        handles.tempRPM = RPMSelected*ones(noPt,1);
        handles.tempEff = yInt;
        handles.tempPR = xNode';
end;
guidata(hObject,handles);


function noPtTurbEdit_Callback(hObject, eventdata, handles)
% hObject    handle to noPtTurbEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noPtTurbEdit as text
%        str2double(get(hObject,'String')) returns contents of noPtTurbEdit as a double


% --- Executes during object creation, after setting all properties.
function noPtTurbEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noPtTurbEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function compInterpMethodBG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compInterpMethodBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in compInterpMethodBG.
function compInterpMethodBG_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in compInterpMethodBG 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function compLinearRB_Callback(hObject, eventdata, handles)
% hObject    handle to compLinearRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compLinearRB


% --- Executes during object creation, after setting all properties.
function radiobutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compPCHPRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function compPCHPRB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compLinearRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compLinearRB


% --- Executes when selected object is changed in compDataBG.
function compDataBG_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in compDataBG 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;

IdxSelected = get(handles.RPMComplistbox,'Value');
RPMSelected = handles.RPMComp_rep{IdxSelected};
set(handles.noPtCompEdit,'String',length(handles.compRPMIdx{IdxSelected}));
SB = get(get(handles.compDataBG,'SelectedObject'),'tag');
cla(handles.CompCondInterp);
switch SB
    case 'compMassFlowRB'
        plot(handles.CompCondInterp,...
            handles.SAE_map_comp.PR(handles.compRPMIdx{IdxSelected}),...
            handles.SAE_map_comp.m_dot(handles.compRPMIdx{IdxSelected})...
            ,'*r');    
    case 'compEfficiencyRB'
        plot(handles.CompCondInterp,...
            handles.SAE_map_comp.PR(handles.compRPMIdx{IdxSelected}),...
            handles.SAE_map_comp.eff(handles.compRPMIdx{IdxSelected})...
            ,'*r');    
end;
guidata(hObject,handles);

% --- Executes when selected object is changed in turbDataBG.
function turbDataBG_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in turbDataBG 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;

IdxSelected = get(handles.RPMTurblistbox,'Value');
RPMSelected = handles.RPMTurb_rep{IdxSelected};
set(handles.noPtCompEdit,'String',length(handles.turbRPMIdx{IdxSelected}));
SB = get(get(handles.turbDataBG,'SelectedObject'),'tag');
cla(handles.TurbCondInterp);
switch SB
    case 'turbMassFlowRB'
        plot(handles.TurbCondInterp,...
            handles.SAE_map_turb.PR(handles.turbRPMIdx{IdxSelected}),...
            handles.SAE_map_turb.m_dot(handles.turbRPMIdx{IdxSelected})...
            ,'*r');    
    case 'turbEffRB'
        plot(handles.TurbCondInterp,...
            handles.SAE_map_turb.PR(handles.turbRPMIdx{IdxSelected}),...
            handles.SAE_map_turb.eff(handles.turbRPMIdx{IdxSelected})...
            ,'*r');    
end;
guidata(hObject,handles);


% --- Executes on button press in TurbCondSavePB1.
function TurbCondSavePB1_Callback(hObject, eventdata, handles)
% hObject    handle to TurbCondSavePB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RPMTemp = handles.SAE_map_turb.RPM;
mDotTemp = handles.SAE_map_turb.m_dot;
PRTemp = handles.SAE_map_turb.PR;
EffTemp = handles.SAE_map_turb.eff;

idx = (RPMTemp == handles.tempTurbRPM(1));
idx0 = find(idx(1:end-1) ~= idx(2:end));
idx1 = idx0(1)+1;
if (length(idx0) == 1)
    if (idx(1) == 1)
        mDotTemp = [handles.tempTurbMDot;mDotTemp(idx1:end)];
        PRTemp = [handles.tempTurbPR;PRTemp(idx1:end)];
        EffTemp = [handles.tempTurbEff;EffTemp(idx1:end)];
        RPMTemp = [handles.tempTurbRPM;RPMTemp(idx1:end)];
    else
        mDotTemp = [mDotTemp(1:idx1-1);handles.tempTurbMDot];
        PRTemp = [PRTemp(1:idx1-1);handles.tempTurbPR];
        EffTemp = [EffTemp(1:idx1-1);handles.tempTurbEff];
        RPMTemp = [RPMTemp(1:idx1-1);handles.tempTurbRPM];
    end;
else
    idx2 = idx0(2);
    mDotTemp = [mDotTemp(1:idx1-1);handles.tempTurbMDot;mDotTemp(idx2+1:end)];
    PRTemp = [PRTemp(1:idx1-1);handles.tempTurbPR;PRTemp(idx2+1:end)];
    EffTemp = [EffTemp(1:idx1-1);handles.tempTurbEff;EffTemp(idx2+1:end)];
    RPMTemp = [RPMTemp(1:idx1-1);handles.tempTurbRPM;RPMTemp(idx2+1:end)];
end;

handles.SAE_map_turb.m_dot = mDotTemp;
handles.SAE_map_turb.PR = PRTemp;
handles.SAE_map_turb.eff = EffTemp;
handles.SAE_map_turb.RPM = RPMTemp;

guidata(hObject,handles);

function TurbCondSavePB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TurbCondSavePB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function compNewNoPt_Callback(hObject, eventdata, handles)
% hObject    handle to compNewNoPt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of compNewNoPt as text
%        str2double(get(hObject,'String')) returns contents of compNewNoPt as a double


% --- Executes during object creation, after setting all properties.
function compNewNoPt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compNewNoPt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in compRedistPB.
function compRedistPB_Callback(hObject, eventdata, handles)
% hObject    handle to compRedistPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newNoPt = (get(handles.compNewNoPt,'String'));
newNoPt = str2double(newNoPt{1});
RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
newRPM = [];
newPR = [];
newMDot = [];
newEff = [];
cla(handles.CompCondInterp);
cla(handles.TurbCondInterp);
for i = 1:length(RPMComp_rep);
    idx = find(handles.SAE_map_comp.RPM == RPMComp_rep{i});
    tempRPM = handles.SAE_map_comp.RPM(idx);
    tempPR = handles.SAE_map_comp.PR(idx);
    tempMDot = handles.SAE_map_comp.m_dot(idx);
    tempEff = handles.SAE_map_comp.eff(idx);
    oldNoPt = length(tempRPM);
    rangePR = max(tempPR) - min(tempPR);
    rangeMDot = max(tempMDot) - min(tempMDot);
    rangeEff = max(tempEff) - min(tempEff);
    lArcPRMDot = sqrt(((tempPR(2:end)-tempPR(1:end-1))/rangePR).^2 + ...
        ((tempMDot(2:end)-tempMDot(1:end-1))/rangeMDot).^2);
    lArcPREff = sqrt(((tempPR(2:end)-tempPR(1:end-1))/rangePR).^2 + ...
        ((tempEff(2:end)-tempEff(1:end-1))/rangeEff).^2);
    totLArcPRMDot = sum(lArcPRMDot);
    totLArcPREff = sum(lArcPREff);
    newPRTemp = zeros(newNoPt,1);
    newMDotTemp = newPRTemp;
    newEffTemp = newPRTemp;
    newRPMTemp = newPRTemp;
    
    newPRTemp(1) = tempPR(1);
    newMDotTemp(1) = tempMDot(1);
    newEffTemp(1) = tempEff(1);
    newRPMTemp(1) = tempRPM(1);
    newPRTemp(newNoPt) = tempPR(oldNoPt);
    newMDotTemp(newNoPt) = tempMDot(oldNoPt);
    newEffTemp(newNoPt) = tempEff(oldNoPt);
    newRPMTemp(newNoPt) = tempRPM(oldNoPt);
    for iii=1:newNoPt-2
        k = 1;
        lArc = totLArcPRMDot/(newNoPt-1)*iii;
        while(lArc>sum(lArcPRMDot(1:k)))
            k = k + 1;
        end;
        newPRTemp(iii+1) = tempPR(k) + (tempPR(k+1) - tempPR(k))* ...
            (lArc - sum(lArcPRMDot(1:k-1)))/lArcPRMDot(k);
    end;
    pp1 = csaps(tempPR,tempMDot);
    pp2 = csaps(tempPR,tempEff);
    newMDotTemp = ppval(pp1,newPRTemp);
    newEffTemp = ppval(pp2,newPRTemp);
    newRPMTemp = tempRPM(1)*ones(newNoPt,1);
    newPR = [newPR; newPRTemp];
    newMDot = [newMDot; newMDotTemp];
    newEff = [newEff; newEffTemp];
    newRPM = [newRPM;newRPMTemp];
    plot(handles.CompCondInterp,tempMDot,tempPR,'o');
    hold(handles.CompCondInterp,'on');
    plot(handles.CompCondInterp,newMDotTemp,newPRTemp,'*');
    plot(handles.TurbCondInterp,tempPR,tempEff,'o');
    hold(handles.TurbCondInterp,'on');
    plot(handles.TurbCondInterp,newPRTemp,newEffTemp,'*');
end;
handles.SAE_map_comp.m_dot = newMDot;
handles.SAE_map_comp.RPM = newRPM;
handles.SAE_map_comp.PR = newPR;
handles.SAE_map_comp.eff = newEff;
guidata(hObject,handles);



function SAE_map_comp_new = produceNewSAEComp(RPMNew,SAE_map_comp)
mDot = SAE_map_comp.m_dot;
PR = SAE_map_comp.PR;
eff = SAE_map_comp.eff;
RPM = SAE_map_comp.RPM;

RPM_rep = unique(RPM);
m = length(RPM_rep);
ppMDot = cell(m,1);
ppEff = ppMDot;
for i = 1:m;
    idx = RPM == RPM_rep(i);
    mDotTemp = mDot(idx);
    PRTemp = PR(idx);
    effTemp = eff(idx);
    RPMTemp = RPM(idx);
    ppMDot = csaps(PRTemp,mDotTemp);
    ppEff = csaps(PRTemp,effTemp);
    PRGrid(i,:) = ppMDot.breaks;
    ppMDotGrid(i,:,:) = ppMDot.coefs;
    ppEffGrid(i,:,:) = ppEff.coefs;
    noOrd = ppMDot.order;
end;
% points interpolation
n = length(PRTemp);
for i = 1:n
    ppPR = csaps(RPM_rep,PRGrid(:,i));
    PRNew(:,i) = ppval(ppPR,RPMNew);
    if i < n
        for ii = 1:noOrd
            ppMDot = csaps(RPM_rep,ppMDotGrid(:,i,ii));
            ppEff = csaps(RPM_rep,ppEffGrid(:,i,ii));
            ppMDotGridNew(:,i,ii) = ppval(ppMDot,RPMNew);
            ppEffGridNew(:,i,ii) = ppval(ppEff,RPMNew);
        end;
    end;
end;
RPMNew = RPMNew*ones(1,n);
n = size(RPMNew);
fig1 = figure(1)
h1 = axes
fig2 = figure(2)
h2 = axes
for i = 1:n
    ppMDot.breaks = PRNew(i,:);
    ppMDot.coefs = permute(ppMDotGridNew(i,:,:),[2 3 1]);
    ppMDot.pieces = length(ppMDot.breaks) - 1;
    mDotNew(i,:) = ppval(ppMDot,PRNew(i,:));
    ppEff.breaks = PRNew(i,:);
    ppEff.coefs = permute(ppEffGridNew(i,:,:),[2 3 1]);
    ppEff.pieces = length(ppEff.breaks) - 1;
    effNew(i,:) = ppval(ppEff,PRNew(i,:));
    plot(h1,PRNew(i,:),mDotNew(i,:),'o');
    hold(h1,'on');
    plot(h2,PRNew(i,:),effNew(i,:),'o');
    hold(h2,'on');
end;
SAE_map_comp_new = SAE_map_comp;
SAE_map_comp_new.RPM = reshape(RPMNew',[],1);
SAE_map_comp_new.m_dot = reshape(mDotNew',[],1);
SAE_map_comp_new.eff = reshape(effNew',[],1);
SAE_map_comp_new.PR = reshape(PRNew',[],1);










% --- Executes on button press in regenSAEPB.
function regenSAEPB_Callback(hObject, eventdata, handles)
% hObject    handle to regenSAEPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RPMTemp = handles.SAE_map_comp.RPM;
RPMNew = [RPMTemp(1); (ceil(RPMTemp(1)/1000)*1000:1000:RPMTemp(end))'; RPMTemp(end)];
handles.SAE_map_comp = produceNewSAEComp(RPMNew,handles.SAE_map_comp);

cla(handles.CompMap);
cla(handles.TurbMap);
line(handles.SAE_map_comp.m_dot,handles.SAE_map_comp.PR,handles.SAE_map_comp.eff,...
    'marker','.','markersize',4,'linestyle','none','parent',...
    handles.CompMap);
line(handles.SAE_map_turb.PR,handles.SAE_map_turb.m_dot,handles.SAE_map_turb.eff,...
    'marker','.','markersize',4,'linestyle','none','parent',...
    handles.TurbMap);
handles.RPMComp_rep = num2cell(unique(handles.SAE_map_comp.RPM));
handles.RPMTurb_rep = num2cell(unique(handles.SAE_map_turb.RPM));
for i = 1:length(handles.RPMComp_rep);
    handles.compRPMIdx{i} = find(handles.SAE_map_comp.RPM == ...
        handles.RPMComp_rep{i});
end;
for i = 1:length(handles.RPMTurb_rep);
    handles.turbRPMIdx{i} = find(handles.SAE_map_turb.RPM == ...
        handles.RPMTurb_rep{i});
end;
set(handles.RPMComplistbox,'String',handles.RPMComp_rep);
set(handles.noPtCompEdit,'String',length(handles.compRPMIdx{1}));
set(handles.RPMTurblistbox,'String',handles.RPMTurb_rep);
set(handles.noPtTurbEdit,'String',length(handles.turbRPMIdx{1}));
cla(handles.CompCondInterp);
cla(handles.TurbCondInterp);
plot(handles.CompCondInterp,handles.SAE_map_comp.PR(handles.compRPMIdx{1}),...
    handles.SAE_map_comp.m_dot(handles.compRPMIdx{1}),'*r');
plot(handles.TurbCondInterp,handles.SAE_map_turb.PR(handles.turbRPMIdx{1}),...
    handles.SAE_map_turb.m_dot(handles.turbRPMIdx{1}),'*r');


guidata(hObject,handles);


% --- Executes on button press in expMap.
function expMap_Callback(hObject, eventdata, handles)
% hObject    handle to expMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','SAE_map_comp',handles.SAE_map_comp);
