%% cmap_eff.m
%
%   Script for polynomial fit of the compressor efficiency. Uses curve fitting toolbox. Saves the
%   coefficients to 'ceff.mat'. Plots the results.
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 05.06.2012


load('SAE_map_comp_ABB_TPL_B.mat')
load('SAE_map_turb_ABB_TPL_B.mat')

n = 100;        % Resolution of the meshgrid
Qmin = 3;       % Min volume flow, for plotting
Qmax = 16;      % Max volume flow, for plotting
PRmin = 1;      % Min pressure ratio, for plotting
PRmax = 5;    % Max pressure ratio, for plotting

eQ = SAE_map_comp.m_dot;
eP = SAE_map_comp.PR;
e = SAE_map_comp.eff;
%%
xlin = linspace(min(eQ),max(eQ),n); 
ylin = linspace(min(eP),max(eP),n);
% xlin = linspace(Qmin,Qmax,n); 
% ylin = linspace(PRmin,PRmax,n);

[X,Y] = meshgrid(xlin,ylin);
%Z = griddata(x,y,z,X,Y,'cubic');

eta = surfaceFitEff(eQ,eP,e);
Z = eta(X,Y);
p = coeffvalues(eta);
save('ceff_ABB_TPL_B.mat', 'p')

%% Plotting


% Plotting contours
figure
%contourcmap('copper',0.2)
v = [0.73 0.77 0.8 0.82 0.83];
[C,h] = contour(X,Y,Z,v,'LineWidth',1.5);
clabel(C,h, 'LabelSpacing', 10000, 'Rotation',0,'FontSize',5)
hold on
scatter(eQ,eP, 'ok')
xlabel('Q_c [m^3/s]')
ylabel('\Pi_c [-]')
grid on
%clabel(C,h, 'manual','FontSize',5) 
%axis([Qmin Qmax PRmin PRmax])

% Plotting surface
figure
surf(X,Y,Z, 'EdgeColor','none')
hold on
%plot3(x,y,z,'.b','MarkerSize',15)
shading interp
zlim([0.68 max(Z(:))])
colormap('hot')
colorbar
caxis([0.68 max(Z(:))])
xlabel('Q_c [m^3/s]')
ylabel('\Pi_c [-]')
zlabel('\eta_c [-]')
hold on
plot3(eQ,eP,e,'.b','MarkerSize',15)

% Surface with non black around when seing two axes
figure
surf(X,Y,Z, 'EdgeColor','none')
hold on
%plot3(x,y,z,'.b','MarkerSize',15)
shading interp
axis([Qmin Qmax PRmin PRmax 0.68 max(Z(:))])
colormap('hot')
colorbar
caxis([0.68 max(Z(:))])
xlabel('Q_c [m^3/s]')
ylabel('\Pi_c [-]')
zlabel('\eta_c [-]')
% hold on
% plot3(eQ,eP,e,'.b','MarkerSize',15)


% % Plotting mesh with samples
% figure
% mesh(X,Y,Z);
% axis tight;
% hold on
% plot3(eQ,eP,e,'.b','MarkerSize',15)
% colormap('hot')
% colorbar
% zlim([0.68 max(Z(:))])
% caxis([0.68 max(Z(:))])
% xlabel('Q_c [m^3/s]')
% ylabel('\Pi_c [-]')
% zlabel('\eta_c [-]')
