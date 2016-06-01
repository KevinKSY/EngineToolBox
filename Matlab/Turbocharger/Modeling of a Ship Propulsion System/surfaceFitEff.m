function [fitresult, gof] = surfaceFitEff(eQ, eP, e)
%SURFACEFITEFF(EQ,EP,E)
%  Fit surface to compressor efficiency data using a 5th order polynomial of both x and y
%
%      X Input : eQ - vector of volume flows [m^3/s]
%      Y Input : eP - vector of pressure ratios [-]
%      Z Output: e - vector of efficiencies [-]
%  Output:
%      fitresult : an sfit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, SFIT.
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 03/2012
%% Fit: 'efficiencyFit'.
[xInput, yInput, zOutput] = prepareSurfaceData( eQ, eP, e );

% Set up fittype and options.
ft = fittype( 'poly55' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( [xInput, yInput], zOutput, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'efficiencyFit' );
% h = plot( fitresult, [xInput, yInput], zOutput );
% legend( h, 'efficiencyFit', 'e vs. eQ, eP', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'eQ' );
% ylabel( 'eP' );
% zlabel( 'e' );
% grid on
% view( -13.5, 32 );
% zlim([0.70 0.85]);
% colormap('jet');
% colorbar
