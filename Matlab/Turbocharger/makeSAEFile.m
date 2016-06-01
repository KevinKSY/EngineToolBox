function makeSAEFile()
filename = 'C:\Users\yum\Documents\GitHub\EngineToolBox\Matlab\Turbocharger\TurbinkartScania.SAE';
[NTurb,dmTurb,ERTurb,etaTurb] = importfile1(filename);
filename = 'C:\Users\yum\Documents\GitHub\EngineToolBox\Matlab\Turbocharger\CompressorMapScania.SAE';
[NComp,dmComp,PRComp,etaComp] = importfile(filename);

nPR = 1000.

SAE_map_comp.model = 'SCANIA_DI09_072M';
SAE_map_comp.inlet_diam = [];
SAE_map_comp.outlet_diam = [];
SAE_map_comp.type = 'radial';
SAE_map_comp.inertia = [];

RPMUniq = unique(NComp);
SAE_map_comp.RPM_surge = zeros(length(RPMUniq),1);
SAE_map_comp.PR_surge = zeros(length(RPMUniq),1);
SAE_map_comp.m_dot_surge = zeros(length(RPMUniq),1);
SAE_map_comp.RPM_chock = zeros(length(RPMUniq),1);
SAE_map_comp.PR_chock = zeros(length(RPMUniq),1);
SAE_map_comp.m_dot_chock = zeros(length(RPMUniq),1);
fig3= figure( 'Name', 'Turbine m_dot vs. PR' );
ax3 = axes('Parent',fig3);
hold(ax3,'on');
fig4 = figure('Name', 'Turbine eff vs. PR' );
ax4 = axes('Parent',fig4);
hold(ax4,'on');

for i = 1:length(RPMUniq)
    idx = find(NComp == RPMUniq(i));
    SAE_map_comp.RPM_surge(i) = RPMUniq(i);
    SAE_map_comp.PR_surge(i) = PRComp(idx(1));
    SAE_map_comp.m_dot_surge(i) = dmComp(idx(1));
    SAE_map_comp.RPM_chock(i) = RPMUniq(i);
    SAE_map_comp.PR_chock(i) = PRComp(idx(end));
    SAE_map_comp.m_dot_chock(i) = dmComp(idx(end));
    PRTemp = PRComp(idx);
    dmTemp = dmComp(idx);
    etaTemp = etaComp(idx);
    ft = 'pchipinterp';
    % Fit model to data.
    [fitresult] = fit( PRTemp, dmTemp, ft, 'Normalize', 'on' );
    dPR = (PRTemp(1) - PRTemp(end))/nPR;
    PRNew = (PRTemp(end):dPR:PRTemp(1))';
    dmNew = fitresult(PRNew);
    
    [fitresult] = fit( PRTemp, etaTemp, ft, 'Normalize', 'on' );
    etaNew = fitresult(PRNew);
    
    plot(ax3,PRNew,dmNew,':');
    plot(ax3,PRTemp,dmTemp,'r*');
    plot(ax4,PRNew,etaNew);
    plot(ax4,PRTemp,etaTemp);
    if i == 1
        SAE_map_comp.RPM = RPMUniq(i)*ones(length(PRNew),1);
        SAE_map_comp.m_dot = dmNew;
        SAE_map_comp.eff = etaNew;
        SAE_map_comp.PR = PRNew;
    else
        SAE_map_comp.RPM = [SAE_map_comp.RPM; RPMUniq(i)*ones(length(PRNew),1)];
        SAE_map_comp.m_dot = [SAE_map_comp.m_dot; dmNew];
        SAE_map_comp.eff = [SAE_map_comp.eff; etaNew];
        SAE_map_comp.PR = [SAE_map_comp.PR; PRNew];
    end;    
end;


dER = 0.05;
SAE_map_turb.model = 'SCANIA_DI09_072M';
SAE_map_turb.type = 'radial'
SAE_map_turb.inertia = [];
SAE_map_turb.eta_max = max(etaTurb);
RPMUniq = unique(NTurb);
fig = figure( 'Name', 'Turbine m_dot vs. PR' );
ax = axes('Parent',fig);
hold(ax,'on');
fig2 = figure('Name', 'Turbine eff vs. PR' );
ax2 = axes('Parent',fig2);
hold(ax2,'on');
for i = 1:length(RPMUniq)
    idx = find(NTurb == RPMUniq(i));
    ERTemp = ERTurb(idx);
    dmTurbTemp = dmTurb(idx);
    etaTurbTemp = etaTurb(idx);
    ft = fittype( 'power2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.024709030737215 0.506599521116278 -0.000390702358136355];
    fit_m_dot_PR = fit( ERTemp, dmTurbTemp, ft, opts );
    % Plot fit with data.
    %axes(ax);
    %h = plot(fit_m_dot_PR,ERTemp,dmTurbTemp);
    err = 1;
    x = ERTemp(1);
    while err > 1e-4
        y = fit_m_dot_PR(x);
        dy = (fit_m_dot_PR(x+x/1000) - fit_m_dot_PR(x-x/1000))/(x/500);
        dx = - y/dy;
        err = abs(dx/x);
        if x + dx < 0 
            x = 0.01*x;
        else
            x = x + dx;
        end;
    end;
    ERTempNew = (x:dER:max(ERTemp))';
    dmTurbTempNew = fit_m_dot_PR(ERTempNew);
    etaTurbTempNew = zeros(length(ERTempNew),1);
    
    plot(ax, ERTempNew,dmTurbTempNew);
    plot(ax, ERTemp,dmTurbTemp,'*');
    for j = 1:length(ERTempNew)
        if ERTempNew(j) < min(ERTemp)
            etaTurbTempNew(j) = (etaTurbTemp(1));
        else
            etaTurbTempNew(j) = interp1(ERTemp,etaTurbTemp,ERTempNew(j),'PCHIP');
        end;
    end;
    plot(ax2,ERTemp,etaTurbTemp,'*');
    plot(ax2,ERTempNew,etaTurbTempNew);
    if i == 1
        SAE_map_turb.RPM = RPMUniq(i)*ones(length(ERTempNew),1);
        SAE_map_turb.m_dot = dmTurbTempNew;
        SAE_map_turb.eff = etaTurbTempNew;
        SAE_map_turb.PR = ERTempNew;
    else
        SAE_map_turb.RPM = [SAE_map_turb.RPM; RPMUniq(i)*ones(length(ERTempNew),1)];
        SAE_map_turb.m_dot = [SAE_map_turb.m_dot; dmTurbTempNew];
        SAE_map_turb.eff = [SAE_map_turb.eff; etaTurbTempNew];
        SAE_map_turb.PR = [SAE_map_turb.PR; ERTempNew];
    end;
end;
    
T_ref = [];
D_wheel = [];


save(['SAE_map_comp_' SAE_map_comp.model '.mat'],'SAE_map_comp');
save(['SAE_map_turb_' SAE_map_turb.model '.mat'],'SAE_map_turb');

TC = turbocharger(SAE_map_comp,SAE_map_turb);
close all

figure
ax = axes;
DrawComp(TC,ax,[],[],1e5,300,[]);
figure
DrawTurb(TC,[],[]);
TC.WriteTCMap();

end

function [NTurb,dmTurb,ERTurb,etaTurb] = importfile1(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [NTURB,DMTURB,ERTURB,ETATURB] = IMPORTFILE(FILENAME) Reads data from
%   text file FILENAME for the default selection.
%
%   [NTURB,DMTURB,ERTURB,ETATURB] = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [NTurb,dmTurb,ERTurb,etaTurb] = importfile('TurbinkartScania.SAE',8, 49);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/04/30 13:49:35

%% Initialize variables.
if nargin<=2
    startRow = 8;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%13f%16f%16f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
NTurb = dataArray{:, 1};
dmTurb = dataArray{:, 2};
ERTurb = dataArray{:, 3};
etaTurb = dataArray{:, 4};
end

function [NComp,dmComp,PRComp,etaComp] = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [EXTRAPOLAT,EDMAPFORLOWF,LOWSBYSIMILARI,TYMETHOD] =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   [EXTRAPOLAT,EDMAPFORLOWF,LOWSBYSIMILARI,TYMETHOD] =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [Extrapolat,edmapforlowf,lowsbysimilari,tymethod] = importfile('CompressorMapScania.SAE',6, 59);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/04/30 13:51:58

%% Initialize variables.
if nargin<=2
    startRow = 6;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%13f%16f%16f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
NComp = dataArray{:, 1};
dmComp = dataArray{:, 2};
PRComp = dataArray{:, 3};
etaComp = dataArray{:, 4};


end

