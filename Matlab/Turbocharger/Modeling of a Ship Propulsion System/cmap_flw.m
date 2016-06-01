%% cmap_flw.m
%
%   Script for modeling interpolation and curve fitting of the compressor flow
%   rate (speed curves).
%
%   Interpolation of speed samples, surface fit of turbocharger speed, speed curve modeling with least
%   squares fit and plotting the compressor map speed curves.
%   -Plots the speed curve samples
%   -Plots interpolation ('pchip') of the speed curve samples
%   -Performs a surface fit (interpolation) of speed curves using curve fitting toolbox
%   -Plots the 3-dimensional surface fit
%   -Performes a nonlinear least squares fit of the data to the model, saves the coefficients to 'flw.mat' 
%   -Plots the modeled speed curves
%   -Performes a nonlinear least squares fit of the data to the alternative model
%   -Plots the alternative modeled speed curves
%   -Plots a polyfit of the dimensionless choke/surge line for modeling of the speed curves
%   -Plots the dimensionless speed curves
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Last revised: 22.05.2012

run init; 
load('cmap.mat');       % imports 'c', 'U2', 'U2str' - for explanation see 'cmapImport.m'

% Limits of the compressor map (plotting regions)
Qmin = 3;               % Min volume flow
Qmax = 16;              % Max volume flow
PRmin = 1;              % Min pressure ratio
PRmax = 6.5;            % Max pressure ratio
res = 0.01;             % Resolution of the interpolation
  
n = size(c,3);          % Number of speed curves
m = size(c,1);          % Number of samples per speed curve
l = length(eta_double); % Number of efficiency curves * 2

%% Plotting samples from the real map

Q = zeros(m, n);        % Corrected compressor volume flow Q = Q(#sample, corr.speed)
PR = zeros(m, n);       % Pressure ratio over the compressor PR = PR(#sample, corr.speed)
Q(:,:) = c(:,1,:);
PR(:,:) = c(:,2,:);        

% Scatters the speedline samples
markertype = ['+','o','*','x','s','d','^','v','p','>','<','h','.'];
figure
for i = 1:n
    hold on
    scatter(Q(:,i), PR(:,i), markertype(i),'k', 'LineWidth', 1);
end
axis([Qmin Qmax PRmin PRmax])
xlabel('Corrected volume flow rate [m^3/s]');
ylabel('Pressure ratio [-]');

%% Making speed curves by interpolating the samples (Piecewise cubic Hermite interpolation)

PRi = NaN((PRmax-PRmin)/res+1,n);
Qi = PRi;
for i = 1:n
    pri = linspace(PR(1,i),PR(m,i),round((PR(1,i)-PR(m,i))/res+1));
    qi = interp1(PR(:,i), Q(:,i), pri, 'pchip');
    PRi(round((PRmax-pri(1))/res):round((PRmax-pri(end))/res),i) = pri';
    Qi(round((PRmax-pri(1))/res):round((PRmax-pri(end))/res),i) = qi';
    plothandle(1) = plot(qi, pri, 'k', 'LineWidth',1.5);
    hold on
    text(Q(1,i),PR(1,i),U2str(i), 'FontSize', 5,'HorizontalAlignment','center','VerticalAlignment','bottom')
end

% Plotting the surge and choke lines
hold on
plothandle(2) = plot(Q(1,:),PR(1,:),'k--');
plothandle(3) = plot(Q(end,:),PR(end,:), 'k-.');
legend(plothandle(1:3), {'Speed curve', 'Surge Line', 'Choke Line'})
title('Interpolated constant speed curves + surge and choke lines')

%% Creating a surface fit
% Is commented out because it needs the curvefitting toolbox. 
% (Making a version independent of the toolbox should be easy if the toolbox is not available)
 
% % Makes a surface fit (speed is z-axis) by sftool. Plots the surface i
% % 3-dimensions
% surfacefit = surfaceFitU(Qvec,PRvec,U2vec);
% 
% % Plotting projection on the (Q,PR)-plane
% hold on
% scatter(Qvec, PRvec, 'kx')
% for i=1:n
%     plot(Qi, PRi, 'r');
% end
% 
% % Plotting the surge and choke lines
% hold on
% plot(Q(1,:),PR(1,:),'k','LineWidth',1);
% plot(Q(end,:),PR(end,:), 'k','LineWidth',1);
% 
% % Surface fit with Q as output - ex. Qsurface(3, 450)
% Qsurface = surfaceFitQ(PRvec, U2vec, Qvec);

%% Fitting data for the speed line model

%Making vectors of corresponding blade speed, volume flow rate and pressure ratio
U2vec = ones(m*n,1);    % Compressor blade speed [m/s]
Qvec = ones(m*n,1);     % Corrected volume flow [m^3/s]
PRvec = ones(m*n,1);    % Pressure ratio [-]
for i=1:n
    Qvec(1+(i-1)*m:i*m) = c(:,1,i);
    PRvec(1+(i-1)*m:i*m) = c(:,2,i);
    U2vec(1+(i-1)*m:i*m) = U2(i);
end

% Performing a nonlinear least squares fit of the data to the model in
% 'compflw'
options = optimset('MaxFunEvals',4000,'MaxIter',1000); 
k0 = ones(15,1);    % Initial coefficients
%k0 = [0 -1 2 -2 1 -59 147 -67 1 -6 13 -13 6 1];
k_flw = lsqcurvefit(@compflw,k0,PRvec,Qvec,[],[],options,U2vec); 
%save('cflw.mat','k_flw'); % Saving the coefficients for later use

%% Plotting the modeled speed lines

% Scatters the speedline samples
markertype = ['+','o','*','x','s','d','^','v','p','>','<','h','.'];
figure
hold on
for i = 1:n
    hold on
    scatter(Q(:,i), PR(:,i), markertype(i),'k', 'LineWidth', 1);
end
axis([Qmin Qmax PRmin PRmax])
xlabel('Corrected volume flow rate [m^3/s]');
ylabel('Pressure ratio [-]');

% Plotting modeled speed curves for all samples
for i=1:n
    PRt = PRvec(i*m):(PRvec(1+(i-1)*m)-PRvec(i*m))/30:PRvec(1+(i-1)*m);
    hold on
    PRt = PRt';
    Qt = compflw(k_flw,PRt,U2(i));
    exp_b = plot(Qt, PRt, 'k');
    text(Q(1,i),PR(1,i),U2str(i), 'FontSize', 5,'HorizontalAlignment','center','VerticalAlignment','bottom')
end
title('Model fit to compressor flow samples')

% Plotting a modeled speed curve with U_c = 425 over all pressure ratios
hold on
plot(compflw(k_flw,PRmin:0.01:PRmax,425), PRmin:0.01:PRmax, 'k', 'LineWidth', 1.5);
text(Q(1,4),PR(1,4),'U_2 = 425', 'FontSize', 5,'HorizontalAlignment','center','VerticalAlignment','bottom')

%hold on
%plot(compflw(k_flw,PRmin:0.01:PRmax,200), PRmin:0.01:PRmax, 'k', 'LineWidth', 1.5);
%text(Q(1,4),PR(1,4),'U_2 = 200', 'FontSize', 5,'HorizontalAlignment','center','VerticalAlignment','bottom')

% Plotting speed curves samples in (Phi,Psi)-coordinates
Psi = zeros(m-1,n);
Phi = zeros(m-1,n);
M = zeros(n,1);
for j = 1:n
    Psi(:,j) = (c_pa*T_amb*(PR(2:end,j).^((k_a-1)/k_a)-1))./(0.5*U2(j)^2);
    Phi(:,j) = Q(2:end,j) ./ ((pi/4)*d_c^2*U2(j));
    M(j) = U2(j) / sqrt(k_a*R_a*T_amb);
end

map = colormap(jet(9));
figure
set(gcf, 'defaultaxescolororder',map)
plot(Phi,Psi,'x', 'LineWidth',1.5)
legend(U2str, 'Location','Best')
xlabel('\Phi')
ylabel('\Psi')

% Plotting the model fit in dimensionless quantities
Psi = zeros(9, length(0.65:0.01:1.3));
Phi_new = zeros(size(Psi));
set(gcf, 'defaultaxescolororder',map)
for i = 1:9
    Psi(i,:) = 0.65:0.01:1.3;
    PRinput = nthroot(1+(0.5*U2(i)^2*Psi(i,:))/(c_pa*T_amb),(k_a-1)/k_a);
    Phi_new(i,:) = compflw(k_flw,PRinput,U2(i)) /((pi/4)*d_c^2.*U2(i));
end
hold on
plot(Phi_new', Psi',':')
title('Dimensionless samples with model fit')


%% Fitting data to the alternative flow rate model  (Moraal & Kolmanovsky)

options = optimset('MaxFunEvals',5000,'MaxIter',1000);
k0 = [1 1 1 -1 -1 -1 1 1 1]; % Initial coefficient guess
k_flw2 = lsqcurvefit(@compflw_alt,k0,PRvec,Qvec,[],[],options,U2vec); 


%% Plotting the alternative model (Moraal & Kolmanovsy)

% Scatters the speedline samples
markertype = ['+','o','*','x','s','d','^','v','p','>','<','h','.'];
figure
hold on
for i = 1:n
    hold on
    scatter(Q(:,i), PR(:,i), markertype(i),'k', 'LineWidth', 1);
end
axis([Qmin Qmax PRmin PRmax])
xlabel('Corrected volume flow rate [m^3/s]');
ylabel('Pressure ratio [-]');
title('Alternative model - fit to compressor flow samples')

% Plotting modeled speed curves for all samples
for i=1:n
    PRt = PRvec(i*m):(PRvec(1+(i-1)*m)-PRvec(i*m))/30:PRvec(1+(i-1)*m);
    hold on
    PRt = PRt';
    Qt = compflw_alt(k_flw2,PRt,U2(i));
    alt = plot(Qt, PRt, 'k');
    text(Q(1,i),PR(1,i),U2str(i), 'FontSize', 5,'HorizontalAlignment','center','VerticalAlignment','bottom')
end

% Plotting the modeled speed curve for U_c = 425 (skipping the asymptote)
hold on
PRt = PRmin:0.01:3.14;
Qt = compflw_alt(k_flw2,PRt,425);
plot(Qt, PRt, 'k--', 'LineWidth', 1.5); % Plotting the lower left part of the hyperbola
PRt = 3.16:0.01:PRmax;
Qt = compflw_alt(k_flw2,PRt,425);
plot(Qt, PRt, 'k--', 'LineWidth', 1.5); % Plotting the upper right part of the hyperbola

%% Plotting the choke line, surge line and speed curves in dimensionless parameters

M = U2 ./ sqrt(k_a*R_a*T_amb);    % Mach speed number
Mp = min(M):0.01:max(M);        % Plotting values of Mach speed

% Nondimensional surge pressure and flow
Psi_srg = (c_pa*T_amb*(PR(1,:).^((k_a-1)/k_a)-1))./(0.5*U2.^2); 
Phi_srg = Q(1,:) ./ ((pi/4)*d_c^2*U2);

% Nondimensional choke pressure and flow
Phi_chk = Q(end,:) ./ ((pi/4)*d_c^2*U2); % Maximum nondimensional flow, i.e. choked flow
Psi_chk = (c_pa*T_amb*(PR(end,:).^((k_a-1)/k_a)-1))./(0.5*U2.^2);

Phi_plot = Phi_chk; % Phi_chk or Phi_srg
Psi_plot = Psi_chk; % Psi_chk or Psi_srg

% Plotting choke flow samples and polyfit for 2nd...4th order polynomial
p = polyfit(M, Phi_plot, 2);
figure
title('Choke line samples and polynomial fits')
subplot(2,3,1)
scatter(M,Phi_plot,'k','LineWidth',1)
hold on
plot(Mp, polyval(p,Mp),'k')
ylabel('\Phi');
%text(0.55, 0.128,'(a)', 'fontweight', 'bold')

p = polyfit(M, Phi_plot, 3);
subplot(2,3,2)
scatter(M,Phi_plot,'k','LineWidth',1)
hold on
plot(Mp, polyval(p,Mp),'k')
%text(0.55, 0.128,'(b)', 'fontweight', 'bold')

p = polyfit(M, Phi_plot, 4);
subplot(2,3,3)
scatter(M,Phi_plot,'k','LineWidth',1)
hold on
plot(Mp, polyval(p,Mp),'k')
%text(0.55, 0.128,'(c)', 'fontweight', 'bold')

% Plotting choke pressure for 3rd...5th order polynomial
p = polyfit(M, Psi_plot, 3);
subplot(2,3,4)
scatter(M,Psi_plot,'k','LineWidth',1)
hold on
plot(Mp, polyval(p,Mp),'k')
ylabel('\Psi');
xlabel('M');
%text(0.55, 0.82,'(d)', 'fontweight', 'bold')

p = polyfit(M, Psi_plot, 4);
subplot(2,3,5)
scatter(M,Psi_plot,'k','LineWidth',1)
hold on
plot(Mp, polyval(p,Mp),'k')
xlabel('M');
%text(0.55, 0.82,'(e)', 'fontweight', 'bold')

p = polyfit(M, Psi_plot, 5);
subplot(2,3,6)
scatter(M,Psi_plot,'k','LineWidth',1)
hold on
plot(Mp, polyval(p,Mp),'k')
xlabel('M');
ylabel('\Phi');
%text(0.55, 0.82,'(f)', 'fontweight', 'bold')

% Plotting speed curves in (Phi,Psi)-coordinates
Psi = zeros(m-1,n);
Phi = zeros(m-1,n);
M = zeros(n,1);
for j = 1:n
    Psi(:,j) = (c_pa*T_amb*(PR(2:end,j).^((k_a-1)/k_a)-1))./(0.5*U2(j)^2);
    Phi(:,j) = Q(2:end,j) ./ ((pi/4)*d_c^2*U2(j));
    M(j) = U2(j) / sqrt(k_a*R_a*T_amb);
end

map = colormap(jet(9));
figure
set(gcf, 'defaultaxescolororder',map)
plot(Phi,Psi,'x-', 'LineWidth',1.5)
title('Dimensionless speed curves')
legend(U2str, 'Location','Best')
xlabel('\Phi')
ylabel('\Psi')

%% Testing the polynomial fit to see extrapolation capabilities

% Plotting the surge line flow samples
figure
scatter(U2,Q(1,:),'k','LineWidth',1.5)

% Making polynomial fits
p_chk = polyfit(U2,Q(1,:),2);
U2i = 0:0.01:700;
hold on
plot(U2i,polyval(p_chk,U2i),'k')
xlim([0 700])
ylim([0 16])
xlabel('Compressor blade tip speed [m/s]')
ylabel('Volume flow rate [m^3/s]')
title('Surge line volume flow as polynomial fit of compressor blade tip speed')

% Plotting the surge line pressure samples
figure
scatter(U2,PR(1,:),'k','LineWidth',1.5)

% Making polynomial fits
p_srg = polyfit(U2,PR(1,:),2);
U2i = 0:0.01:700;
hold on
plot(U2i,polyval(p_srg,U2i),'k')
xlim([0 700])
ylim([1 6])
xlabel('Compressor blade tip speed [m/s]')
ylabel('Pressure ratio')
title('Surge line pressure ratio as polynomial fit of compressor blade tip speed')

%% Plotting exponential functions for illustrations in the model development

figure
set(gcf,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|-.|--|:')
y = 0:0.01:3; % Pressure ratio

subplot(2,2,1)
plot(11-exp(y),y,9-exp(y),y)
xlim([0 10])
ylim([0 3])
legend('\Phi = 11-e^{\Psi}','\Phi = 9-e^{\Psi}','Location','SouthWest');
xlabel('\Phi')
ylabel('\Psi')
text(9,2.75,'(a)', 'fontweight', 'bold')

subplot(2,2,3)
plot(11-exp(y),y,11-exp(y-0.5),y)
xlim([0 10])
legend('\Phi = 11-e^{\Psi}','\Phi = 11-e^{\Psi-0.5}','Location','SouthWest');
xlabel('\Phi')
ylabel('\Psi')
text(9, 2.75,'(b)', 'fontweight', 'bold')

subplot(2,2,[2 4])
plot(11-exp(y),y,(11-exp(y.^2)),y);
xlim([0 10])
ylim([0 3])
legend('\Phi = 11-e^{\Psi}','\Phi = 11-e^{\Psi^2}','Location','SouthWest');
xlim([0 10])
xlabel('\Phi')
ylabel('\Psi')
text(9, 2.85,'(c)', 'fontweight', 'bold')