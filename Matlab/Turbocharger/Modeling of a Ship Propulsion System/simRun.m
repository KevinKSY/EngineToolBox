%% engineRun.m
%
%   Simulates the engine in simulink and plots the results
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Last revised:  31.05.2012
%%

clear all;
clc;

run init;                           % Initializes the engine variables
load('cflw');                       % Loads 'k_flw', compressor flow model coefficients
load('ceff');                       % Loads 'p', compressor efficiency polynomial coefficients
load('tflw');                       % Loads 'k_trb', turbine flow model coefficients

simtime = 1000;                      % Total simulation time [s]

% Setting simulation parameters
paramNameValStruct.AbsTol         = '1e-10'; 
paramNameValStruct.RelTol         = '1e-10';
%paramNameValStruct.SimulationMode = 'rapid';
%paramNameValStruct.SaveState      = 'on';
%paramNameValStruct.StateSaveName  = 'xoutNew';
%paramNameValStruct.SaveOutput     = 'on';
%paramNameValStruct.OutputSaveName = 'youtNew';

% Running simulation in simulink
simout = sim('simDiagram.mdl', paramNameValStruct); 

%% Importing variables from simulation

t = simout.get('time');             % Simulation time
states = simout.get('states');      % Simulation states [p_im, T_im, p_em, T_em, w_tc, w_e]
interm = simout.get('interm');      % Intermediate / plotting variables
u = simout.get('u');                % Fuel rack
T_p = simout.get('T_p');            % Propeller thrust
Q_p = simout.get('Q_p');            % Propeller torque
w = simout.get('w');                % Wake coefficient  
V_s = simout.get('V_s');            % Ship surge velocity [m/s]
U_PIb = simout.get('U_PIb');        % Unbounded PI fuel index
N_d = simout.get('N_d');            % Desired engine speed [rpm]
T_ext = simout.get('T_ext');        % External resistance on the ship

% Rename simulation variables
p_im = states(:,1);                 % Pressure intake manifold [Pa]
T_im = states(:,2);                 % Temperature intake manifold [K]
p_em = states(:,3);                 % Pressure exhaust manifold [Pa]
T_em = states(:,4);                 % Temperature exhaust manifold [K]
w_tc = states(:,5);                 % Turbocharger angular frequency [rad/s]
w_e = states(:,6);                  % Engine angular frequency [rad/s]

N_e = (60/(2*pi))*w_e;              % Engine rotational speed [rpm]
N_tc = (60/(2*pi))*w_tc;            % Turboshaft rotational speed [rpm]

T_c = interm(:,1);                  % Temperature of air after compressor [k]
mdot_c = interm(:,2);               % Compressor mass flow rate [kg/s]
mdot_a = interm(:,3);               % Air mass flow rate into engine cylinders [kg/s]
mdot_f = interm(:,5);               % Injected fuel mass flow rate [kg/s]
mdot_e = interm(:,6);               % Exhaust mass flow rate from cylinders [kg/s]
AFR = interm(:,7);                  % Air-to-fuel ratio
eta_th = interm(:,8);               % Combustion efficiency
p_i = interm(:,9);                  % Indicated mean effective pressure [Pa]
p_f = interm(:,10);                 % Friction mean effective pressure [Pa]
p_e = interm(:,11);                 % Break mean effective pressure [Pa]
zeta = interm(:,12);                % Fuel chemical energy proportian in exhaust gas
T_e = interm(:,13);                 % Exhaust gas (from cylinders) temperature [K]
p_ep = interm(:,14);                % Pressure in exhaust pipe / turbine outlet [Pa]
T_ep = interm(:,15);                % Temperature in exhaust pipe / turbine outlet [K]
mdot_t = interm(:,16);              % Turbine mass flow rate [kg/s]
Q_c = interm(:,17);                 % Compresser absorbed torque [Nm]
Q_t = interm(:,18);                 % Turbine produced torque [Nm]
Q_e = interm(:,19);                 % Engine produced torque [Nm]
Q_l = interm(:,20);                 % Propeller load torque [Nm]
bsr = interm(:,21);                 % Blade speed ratio
eta_t = interm(:,22);               % Turbine isentropic efficiency
eta_c = interm(:,23);               % Compressor isentropic efficiency
Psi_e = interm(:,24);               % Flow function value for air through engine cylinders
T_i = interm(:,25);                 % Temperature of air after intercooler [K]

% Calculation of ship and propeller variables
R_T = Chat_T * V_s.^2;              % Ship Total Resistance [Nm]
V_knots = (3600/1852)*V_s;          % Convertion from m/s to knots
n = w_e / (2*pi);                   % Propeller speed [rps]
V_a = (1-w).*V_s;                   % Advance speed [m/s]
J_a = V_a ./ (n*D);                 % Advance ratio [-]
eta_0 = (T_p.*V_a)./(2*pi.*n.*Q_p); % Open water efficiency [-]
K_Ti = interp1(J,K_T,J_a,'linear'); % Propeller thrust coefficient [-]
K_Qi = interp1(J,K_Q,J_a,'linear'); % Propeller torque coefficient [-]
S_R = 1 - V_a ./ (PD*D*n);          % Real slip ratio [-]

%% Plotting properties

set(0,'DefaultLineLineWidth',1.5) 
set(0,'DefaultAxesLineWidth',1)
%set(0,'DefaultAxesColorOrder',[0 0 0])
%set(0,'DefaultAxesLineStyleOrder','-|-.|--|:')

%% Plotting engine/turbocharger variables

Uc_low = (250*60)/(pi*d_c);     % Lowest speed of compressor in map
Uc_high = (550*60)/(pi*d_c);    % Highest speed of compressor in map

figure
subplot(5,2,1)
plot(t,u,t,U_PIb,'--')
xlabel('Time [s]')
ylabel('Fuel Rack [-]')
legend('Real','PI','Location','SouthWest')
grid on

subplot(5,2,2)
plot(t, N_e, t, N_mcr*ones(length(t),1), '--r', t, N_d,'-.');
grid on
xlabel('Time [s]')
ylabel('Engine Speed [rpm]')
legend('Engine', 'MCR', 'Setpoint','Location', 'SouthWest')

subplot(5,2,3)
plot(t, p_im/p_amb, t, p_em./p_ep);
grid on
xlabel('Time [s]');
ylabel('Pressure Ratio [-]')
legend('compressor', 'turbine', 'Location', 'SouthWest')

subplot(5,2,4);
plot(t,N_tc, t, ones(length(t),1)*Uc_low, '--r', t,ones(length(t),1)*Uc_high,'--r')
grid on
xlabel('Time [s]')
ylabel('Turbocharger Speed [rpm]')
legend('TC','mapped region', 'Location', 'SouthWest')

subplot(5,2,5)
plot(t, T_c,t, T_im, t, T_em);
grid on
legend('T_c','T_{im}', 'T_{em}', 'Location','SouthWest')
xlabel('Time [s]');
ylabel('Temperature [K]')

subplot(5,2,6)
plot(t, p_e/100000, t, p_i/100000)
grid on
legend('bmep','imep', 'Location', 'SouthWest')
xlabel('Time [s]')
ylabel('Pressure [bar]')

subplot(5,2,7)
plot(t, mdot_c, t, mdot_t)
legend('compressor','turbine', 'Location', 'SouthWest')
xlabel('Time [s]')
ylabel('Mass flow rate [kg/s]')
grid on

subplot(5,2,8)
plot(t, eta_c, t, eta_t)
xlabel('Time [s]')
ylabel('Efficiency [-]')
legend('compressor', 'turbine','Location','SouthWest')
grid on

subplot(5,2,10)
plot(t,AFR)
xlabel('Time [s]')
ylabel('Air-fuel-ratio [-]')
grid on
ylim([35 45])

subplot(5,2,9)
plot(t, Q_e, t, Q_l)
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('Engine','Propeller','Location','SouthWest')
grid on

%% Plotting ship variables

figure
subplot(3,2,1)
plot(t, V_knots)
xlabel('Time [s]')
ylabel('Ship Velocity [knots]')
grid on

subplot(3,2,2)
plot(t,T_p*(1-t_d)/1000,t,R_T/1000 + T_ext/1000,t,T_ext/1000)
xlabel('Time [s]')
ylabel('Thrust, resistance [kN]')
grid on
legend('Thrust to Hull','Total Resistance','Extra Resistance','Location', 'SouthEast')


subplot(3,2,3)
plot(t,J_a)
xlabel('Time [s]')
ylabel('Advance Ratio [-]')
grid on

subplot(3,2,4)
plot(t,eta_0)
xlabel('Time [s]')
ylabel('Open Water Efficiency [-]')
grid on

subplot(3,2,5)
plot(t, K_Ti, t, 10*K_Qi)
xlabel('Time [s]')
ylabel('KT, 10KQ [-]')
grid on
legend('K_T', '10K_Q', 'Location', 'SouthWest')
ylim([0.05 0.2])

subplot(3,2,6)
plot(t, S_R)
xlabel('Time [s]')
ylabel('Real Slip Ratio [-]')
grid on