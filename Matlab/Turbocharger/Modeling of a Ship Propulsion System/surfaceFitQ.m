function [fitresult, gof] = surfaceFitQ(PRvec, U2vec, Qvec)
%SURFACEFITQ(PRVEC,U2VEC,QVEC)
%  Fit surface to speed line samples using cubic interpolation, with Q as output.
%
%  Data for 'Qinterp' fit:
%      X Input : PR - pressure ratio vector [-]
%      Y Output: U_c - blade tip speed vector [m/s]
%      Z Output: Q - volume flow rate vector [m^3/s]
%  Output:
%      fitresult : an sfit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, SFIT.
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 03/2012

%% Fit: 'Qinterp'.
[xInput, yInput, zOutput] = prepareSurfaceData(PRvec,U2vec,Qvec );

% Set up fittype and options.
ft = 'cubicinterp';
opts = fitoptions( ft );
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( [xInput, yInput], zOutput, ft, opts );

end