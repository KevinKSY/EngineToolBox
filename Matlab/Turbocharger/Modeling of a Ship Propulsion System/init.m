%% init.m 
%
%   Script for setting the engine constants and initial values for simulation
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Last revised: 18.06.2012
%%

clear all;
clc;

global k_a R_a T_amb c_pa d_c % for use in 'compflw' and 'compflw_alt', global to avoid passing as parameters

nx = 6;                 % Number of engine states
ni = 25;                % Number of intermediate/plotting variables

% Values for MAN-B&W 6L60MC-C MK6(Xiros 2002, p.34)
P_mcr = 9177*10^3;      % Engine MCR power[W]
N_mcr = 114.6;          % Engine MCR speed[rpm]
z_c = 6;                % Number of engine cylinders
V_d = 0.5429;           % Engine displaced volume pr cylinder [m^3]
J_tot = 59800;          % Engine + propeller nominal moment of inertia [kg m^2]
p_imax = 16*10^5;       % Maximum imep [Pa]
T_ep = 273+230;         % Exhaust pipe temperature [K]
k_f2 = 0.03;                
k_f1 = 681.3110;
k_f0 = 818.2450;
k_f = [k_f2 k_f1 k_f0]; % Friction pressure coefficients

% Thermodynamic values 
T_amb = 293.2;          % Ambient temperature [K]
p_amb = 101325;         % Ambient/athmospheric pressure [Pa]
c_pa = 1005;            % Specific heat capacity air, constant pressure [J/(kg K)]
c_va = 718;             % Specific heat capacity for air at constant volume [J/(kg K)] 
c_pe = 1117;            % Specific heat capacity exhaust gas at constant pressure [J/(kg K)]
k_a = 1.4;              % Ratio of specific heats for air  [J/(kg K]
k_e = 1.34;             % Ratio of specific heats for exhaust [J/(kg K] at (400-1100K)
R_a = 287;              % Specific gas constant air [J/(kg K)] 
R_e = 277.1;            % Specific gas constant exhaust [J/(kg K)]

% Fuel Values
AFR_s = 14.7;           % Stoichiometric air/fuel ratio
AFR_high = 22;          % Air/fuel ratio for perfect combustion,
AFR_low = 8.0;          % Lean flamability limit
H_l = 40640*10^3;       % Fuel lower heating value [J/kg]
m_fmax = 0.0393;        % Max fuel amount per cylinder per cycle [kg] 
k_z1 = 0.0105*10^-5;  
k_z0 = 0.3120;
k_z = [k_z1 k_z0];      % Fuel chemical energy proportion in ex. gas coefficient

% Valves / Orifices
c_d = 0.9;              % Inlet port/exhaust valve flow discharge coefficient
A_Veq = 0.04;           % Inlet port/exhaust valve mean equivalent area [m^2]
Pi_tr = 0.95;           % Treshold for pressure ratio so that calculation retains Lipchitz

% Intercooler
T_water = 298;          % Coolant water temperature [K]
eta_ic = 0.98;          % Intercooler efficiency

% Intake/exhaust manifold
V_im = 1.5*V_d;         % Volume intake manifold [m^3]
V_em = 1.5*V_d;         % Volume exhaust manifold [m^3]
eta_em = 1;             % Exhaust manifold thermodynamic efficiency

% Turbocharger
J_tc = 1.5*4.83;        % Inertia of turbocharger system [kg m^2]
d_c = 0.5;              % Compressor diameter [m]
d_t = 0.5;              % Turbine diameter [m]
eta_t_max = 0.82;       % Turbine max efficiency
bsr_max = 0.7;          % Turbine blade speed ratio at max efficiency
k_t = 2;                % Turbine efficiency coefficient
T_a0 = 293.2;           % Turbine inlet temperature in map [K]
k_srg = [4.5895*10^-5, 0.0018, -0.0645]; % Polynomial coefficients for the surge volume flow

% Propeller load dynamics
rho_a = 1.23;           % Density of air [kg m^-3]
rho_w = 1025;           % Sea water density [kg/m^3]
w = 0.277;              % Wake fraction number
Q_mcr = P_mcr/(((2*pi)/60)*N_mcr);
k_Q0 = Q_mcr / N_mcr^2; % Propeller curve torque coefficient [N m/rpm^2]

% Propeller parameters
D = 6.5;                % Propeller diameter [m]
Z = 5;                  % Number of propeller blades [-]
PD = 0.665;             % Propellerpitch to diameter ratio [-]
EAR = 0.57;             % Propeller area ratio [-]
J_ew = ((0.0703*PD^2*EAR^2)/(pi*Z))*rho_w*D^5; % Entrained water inertia [kg m^2]
J_tot = J_tot + J_ew;   % Shafting system total moment of inertia


% Ship parameters
s_dw = 55000;               % Ship size (at scantligng draught) [dwt]
l_oa = 190;                 % Length overall [m]
l_pp = 183;                 % Length between perpendiculars [m]
B = 32.26;                  % Breadth [m]
D_s = 12.7;                 % Scantling draught [m]
D_d = 11.5;                 % Design draught [m]
M_s = 6.358*10^6;           % Ship mass [kg]
M_a = 0.25*M_s;             % Added mass [kg]
t_d = 0.2;                  % Thrust deduction coefficient [-]

% Governor parameters
Kp = 1;                     % Proportional gain (PID)
Ki = 0.1;                   % Integrator gain (PID)
Kw = 1;                     % Integrator anti-windup gain (inverse of anti-windup reset time constant)
u_lb = 0;                   % Fuel rack lower bound
u_ub = 1;                   % Fuel rack upper bound
k_trq = [0.4, 1, 0.4, 0.8]; % Torqe limiter constants
k_scv = [0.5, 1, 1.2, 2.0]; % Scavenge limiter constants

% Initial values
p_im0 = 3.5*p_amb;          % Intake manifold pressure [Pa]
T_im0 = 300;                % Intake manifold temperature [K]
p_em0 = 3*p_amb;            % Exhaust manifold pressure [Pa]
T_em0 = 600;                % Exhaust manifold temperature [K]
w_e0 = 114.6*(2*pi)/60;     % Engine angular frequency [rad/s]
w_tc0 = 20000*(2*pi)/60;    % Turbocharger angular frequency [rad/s]
u_0 = 1;                    % Initial fuel rack position
V_s0 = 15*0.5144;           % Initial ship speed [m/s]

x0 = [p_im0; T_im0; p_em0; T_em0; w_tc0; w_e0];

% Creating the thrust (K_T) and torque (K_Q) coefficients
J = 0:0.01:0.74;       % Range of advance speed
[C_T,s,t,u,v] = textread('KT.txt', '%f %d %d %d %d'); %Import the polynomial coefficients
K_T = zeros(size(J));
for i = 1:length(K_T)
    for j = 1:length(C_T)
        K_T(i) = K_T(i) + C_T(j)*J(i)^s(j)*PD^t(j)*EAR^u(j)*Z^v(j);
    end
end
[C_Q,s,t,u,v] = textread('KQ.txt', '%f %d %d %d %d'); %Import the polynomial coefficients
K_Q = zeros(size(J));
for i = 1:length(K_Q)
    for j = 1:length(C_Q)
        K_Q(i) = K_Q(i) + C_Q(j)*J(i)^s(j)*PD^t(j)*EAR^u(j)*Z^v(j);
    end
end

% Calulation of ship resistance coefficient
Vmax = 15*0.5144444; 
KTT = ((1-w)*Vmax) / (N_mcr*D/60);
KTTi = interp1(J,K_T,KTT,'linear');
KQQi = interp1(J,K_Q,KTT,'linear');
Chat_T = (1-t_d)*KTTi*rho_w*D^4*(N_mcr/60)^2/Vmax^2;

% Setting default plotting properties
set(0,'defaultlinelinewidth',1.5) 
set(0,'defaultaxeslinewidth',1)
%set(0,'defaultMarkerSize',10)
%set(0,'defaultpatchlinewidth',1)
%set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|-.|--|:')
