%% tmap.m
% 
%   Script that imports turbine flow map values, 
%   performs least squares model fit,
%   and plots examples of turbine flow with varying exhaust manifold temp.
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 05/2012

k_e = 1.34; % Ratio of specific heats for exhaust [J/(kg K] at (400-1100K)

PR_t = [1.26 1.4 1.8 2.2 2.6 3.0 3.4];              % Turbine expansion ratio (p_ep/p_em) 
Q_t = [0.282 0.325 0.367 0.378 0.3825 0.384 0.385]; % Turbine volume flow rate [m^3/(s*sqrt(K))] 
%Q_t = [0.270 0.3125 0.354 0.366 0.370 0.372 0.373]; % Turbine volume flow rate [m^3/(s*sqrt(K))] with smaller nozzle ring
T_em = [600 800 1000];  % Exhaust manifold temperatures for plot
T_emstr = {'T_{em} = 600[K]', 'T_{em} = 800[K]', 'T_{em} = 1000[K]'};

options = []; 
k0 = [0.5; -2.2];   % Initial coefficient guess
k_trb = lsqcurvefit(@turbflw,k0,PR_t,Q_t); 

k0 = ones(2,1);     % Initial coefficient guess
k_trb2 = lsqcurvefit(@(k0,PR_t)turbflw2(k0,PR_t,k_e),k0,PR_t,Q_t);
%save('tflw.mat', 'k_trb', 'k_trb2') % Saving coefficients for later use

% Plotting the samples and the model fit
figure
set(gcf,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|--|-.|:')
scatter(PR_t, Q_t, 'LineWidth', 1.5)
hold on
P = linspace(1.1,3.70,1000);
ph = plot(P, turbflw(k_trb,P),P, turbflw2(k_trb2,P,k_e));
ylabel('Flow rate [m^3/(s\surd{K})]')
xlabel('Pressure ratio [-]')
ylim([0.18 0.4])
legend(ph,'c_3 \surd{1-\Pi_t^{c_4}}','c_0 \Psi','Location','SouthEast')

% Plotting the modeled turbine flow for different turbine inlet
% temperatures
figure
set(gcf,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|--|-.|:')
P = linspace(1,4,100);
plot(P, turbflw(k_trb,P)*sqrt(T_em(1)), P, turbflw(k_trb,P)*sqrt(T_em(2)), P, turbflw(k_trb,P)*sqrt(T_em(3)))
legend(T_emstr, 'Location', 'SouthEast')
ylabel('Volume flow rate [m^3/s]')
xlabel('Pressure ratio [-]')


