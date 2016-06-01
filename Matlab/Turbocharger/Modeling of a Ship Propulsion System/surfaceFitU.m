function [fitresult, gof] = surfaceFitU(Q, PR, U_c)
%SURFACEFITU(Q,PR,U_C)
%  Fit surface to speed line samples using cubic interpolation, with U as output. Plots
%  the results.
%
%  Data:
%      X Input : Q - volume flow rate vector [m^3/s]
%      Y Input : PR - pressure ratio vector [-]
%      Z Output: U_c - blade tip speed vector [m/s]
%  Output:
%      fitresult : an sfit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, SFIT.
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 03/2012

%% Fit: 'cmapfit'.
[xInput, yInput, zOutput] = prepareSurfaceData(Q,PR,U_c);

% Set up fittype and options.
ft = 'cubicinterp';
opts = fitoptions( ft );
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( [xInput, yInput], zOutput, ft, opts );

% Plot fit with data.
figure( 'Name', 'cmapfit' );
h = plot( fitresult, [xInput, yInput], zOutput );
legend( h, 'Map Fit', 'Sample points', 'Location', 'NorthEast' );
% Label axes
xlabel( 'Q [m^3/s]' );
ylabel( 'PR [-]' );
zlabel( 'U_c [m/s]' );
grid on
view( -22, 18 );

ylim([1 5])
end