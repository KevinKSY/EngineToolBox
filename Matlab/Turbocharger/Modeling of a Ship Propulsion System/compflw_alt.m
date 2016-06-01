function Q = compflw_alt(k,PR,U_c)
% Q = COMPFLW_ALT(k,U_c,PR)
%   Calculates the volume flow rate for the alternative volume flow model. 
%
%   Data:
%       k   : vector of coefficients
%       PR  : pressure ratio 
%       U_c : the blade tip speed [m/s]
%   Output:
%       Q   : corrected volume flow [m^3/s]
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 05/2012


global k_a R_a T_amb c_pa d_c

k11 = k(1);
k12 = k(2);
k13 = k(3);
k21 = k(4);
k22 = k(5);
k23 = k(6);
k31 = k(7);
k32 = k(8);
k33 = k(9);

% Coefficients for less accurate version
% k11 = k(1);
% k12 = k(2);
% k21 = k(3);
% k22 = k(4);
% k31 = k(5);
% k32 = k(6);

M = U_c ./ sqrt(k_a*R_a*T_amb);                           % Inlet Mach number
Psi = (c_pa*T_amb*(PR.^((k_a-1)/k_a)-1))./(0.5*U_c.^2);   % Dimensionless head parameter

k1 = k11 + k12.*M + k13.*M.^2;
k2 = k21 + k22.*M + k23.*M.^2;
k3 = k31 + k32.*M + k33.*M.^2;

% Less accurate version
% k1 = k11 + k12.*M;
% k2 = k21 + k22.*M;
% k3 = k31 + k32.*M;

Phi = (k3.*Psi - k1) ./ (k2 + Psi);

Q = Phi*0.25*pi*d_c^2.*U_c;