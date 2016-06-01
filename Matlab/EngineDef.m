function varargout = EngineDef(varargin)
% ENGINEDEF MATLAB code for EngineDef.fig
%      ENGINEDEF, by itself, creates a new ENGINEDEF or raises the existing
%      singleton*.
%
%      H = ENGINEDEF returns the handle to a new ENGINEDEF or the handle to
%      the existing singleton*.
%
%      ENGINEDEF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENGINEDEF.M with the given input arguments.
%
%      ENGINEDEF('Property','Value',...) creates a new ENGINEDEF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EngineDef_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EngineDef_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EngineDef

% Last Modified by GUIDE v2.5 07-Apr-2016 17:52:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EngineDef_OpeningFcn, ...
                   'gui_OutputFcn',  @EngineDef_OutputFcn, ...
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


% --- Executes just before EngineDef is made visible.
function EngineDef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EngineDef (see VARARGIN)
handles.noField = 11;
% Choose default command line output for EngineDef
handles.output = hObject;
handles.uitable1.ColumnName = '';
handles.uitable1.RowName = {'Number of cylinder';'Number of stroke';'Stroke(mm)';'Bore(mm)';'Compression ratio';'a(mm)';'l(mm)';'\phi@EVO';'\phi@EVC';'\phi@IVO';'\phi@IVC'};
% Update handles structure
handles.uitable1.Data = zeros(handles.noField ,1);
handles.uitable1.ColumnEditable = true;
guidata(hObject, handles);


matlabImage = imread('cylDim.png');
image(matlabImage)
axis(handles.axes1,'off')
axis(handles.axes1,'image')
% UIWAIT makes EngineDef wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EngineDef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data = handles.uitable1.Data;
handles.uitable1.RowName = {'Number of cylinder';'Number of stroke';'Stroke(mm)';'Bore(mm)';'Compression ratio';'a(mm)';'l(mm)';'\phi@EVO';'\phi@EVC';'\phi@IVO';'\phi@IVC'};
eng.nStroke = handles.data(2);
eng.nCyl = handles.data(1);
for i = 1:eng.nCyl
    eng.cyl(i).dim.S = handles.data(3)/1000;
    eng.cyl(i).dim.B = handles.data(4)/1000;
    eng.cyl(i).dim.CR = handles.data(5);
    eng.cyl(i).dim.lambda = handles.data(6)/handles.data(7);       %Ratio of length of crank rod to connecting rod 
    eng.cyl(i).dim.vDisp = eng.cyl(i).dim.B^2*pi/4*eng.cyl(i).dim.S;
    eng.cyl(i).dim.vComp = eng.cyl(i).dim.vDisp / (eng.cyl(i).dim.CR - 1);
    eng.cyl(i).exhVVProf.liftUpRef = [0;0.007;0.021;0.046;0.07;0.105;0.161; ...
        0.211;0.267;0.323;0.372;0.435;0.491;0.547;0.604;0.653;0.702;0.744;...
        0.782;0.818;0.853;0.884;0.909;0.93;0.951;0.965;0.979;0.989;0.996; ...
        1];                 %Exhaust normalized lift array for opening [0~1]
    eng.cyl(i).exhVVProf.liftDownRef = [1;0.996;0.993;0.982;0.972;0.958;0.937;...
        0.912;0.884;0.856;0.821;0.786;0.744;0.702;0.656;0.614;0.565;0.509;...
        0.442;0.379;0.319;0.253;0.196;0.14;0.095;0.06;0.032;0.018;0.006;...
        0];                 %Exhaust normalized lift array for closing [0~1]
    eng.cyl(i).exhVVProf.cAStartRef = handles.data(8);   %Nominal crank angle at start of valve lif [deg]
    eng.cyl(i).exhVVProf.dCALiftUp = (handles.data(9) - handles.data(8))/2;     %Nominal duration of crank angle for valve opening [deg]
    eng.cyl(i).exhVVProf.dCALiftDown = (handles.data(9) - handles.data(8))/2;   %Nominal duration of crank angle for valve closing [deg]
    eng.cyl(i).exhVVProf.dCALiftTopRef = 0; %Nominal duration of crank angle for valve open [deg]
    
    eng.cyl(i).intakeVVProf.liftUpRef = [0;0.007;0.021;0.046;0.07;0.105;0.161; ...
        0.211;0.267;0.323;0.372;0.435;0.491;0.547;0.604;0.653;0.702;0.744;...
        0.782;0.818;0.853;0.884;0.909;0.93;0.951;0.965;0.979;0.989;0.996; ...
        1];                    %intake normalized lift array for opening [0~1]
    eng.cyl(i).intakeVVProf.liftDownRef = [1;0.996;0.993;0.982;0.972;0.958;0.937;...
        0.912;0.884;0.856;0.821;0.786;0.744;0.702;0.656;0.614;0.565;0.509;...
        0.442;0.379;0.319;0.253;0.196;0.14;0.095;0.06;0.032;0.018;0.006;...
        0];                 %intake normalized lift array for closing [0~1]
    eng.cyl(i).intakeVVProf.cAStartRef = handles.data(10);   %Nominal crank angle at start of valve lif [deg]
    eng.cyl(i).intakeVVProf.dCALiftUp = (handles.data(11) - handles.data(10))/2;     %Nominal duration of crank angle for valve opening [deg]
    eng.cyl(i).intakeVVProf.dCALiftDown = (handles.data(11) - handles.data(10))/2;   %Nominal duration of crank angle for valve closing [deg]
    eng.cyl(i).intakeVVProf.dCALiftTopRef = 0; %Nominal duration of crank angle for valve open [deg]
    phi = eng.cyl(i).exhVVProf.cAStartRef;
    eng.cyl(i).vEVO = getVolCylRecip(phi,eng.cyl(i).dim.B,eng.cyl(i).dim.S,eng.cyl(i).dim.CR,eng.cyl(i).dim.lambda,true);
    phi = eng.cyl(i).exhVVProf.cAStartRef + eng.cyl(i).exhVVProf.dCALiftUp + eng.cyl(i).exhVVProf.dCALiftDown + eng.cyl(i).exhVVProf.dCALiftTopRef;
    eng.cyl(i).vEVC = getVolCylRecip(phi,eng.cyl(i).dim.B,eng.cyl(i).dim.S,eng.cyl(i).dim.CR,eng.cyl(i).dim.lambda,true);
    phi = eng.cyl(i).intakeVVProf.cAStartRef;
    eng.cyl(i).vIVO = getVolCylRecip(phi,eng.cyl(i).dim.B,eng.cyl(i).dim.S,eng.cyl(i).dim.CR,eng.cyl(i).dim.lambda,true);
    phi = eng.cyl(i).intakeVVProf.cAStartRef + eng.cyl(i).exhVVProf.dCALiftUp + eng.cyl(i).exhVVProf.dCALiftDown + eng.cyl(i).exhVVProf.dCALiftTopRef;
    eng.cyl(i).vIVC = getVolCylRecip(phi,eng.cyl(i).dim.B,eng.cyl(i).dim.S,eng.cyl(i).dim.CR,eng.cyl(i).dim.lambda,true);
end;
assignin('base','eng',eng);