function Q = compflw(k,PR,U_c)
% Q = VOLFLOW(PR,U_c,k)
%   Calculates the volume flow rate for the modeled compressor. 
%
%   Data:
%       PR  : pressure ratio over the compressor (T_im / T_amb) 
%       U_c : blade tip speed of the compressor [m/s]
%       k   : coefficients of the model
%   Output:
%       Q   : air volume flow rate [m^3/s]
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 05/2012

global k_a R_a T_amb c_pa d_c

a4 = k(1);
a3 = k(2);
a2 = k(3);
a1 = k(4);
a0 = k(5);
b2 = k(6);
b1 = k(7);
b0 = k(8);
c5 = k(9);
c4 = k(10);
c3 = k(11);
c2 = k(12);
c1 = k(13);
c0 = k(14);

M = U_c ./ sqrt(k_a*R_a*T_amb);  % Inlet Mach number
Psi = (c_pa*T_amb*(PR.^((k_a-1)/k_a)-1))./(0.5*U_c.^2); % Dimensionless head parameter

% alpha = a4*M.^4 + a3*M.^3 + a2*M.^2 + a1*M + a0;
% beta = b2*M.^2 + b1*M + b0;
% Psi_th = c5*M.^5 + c4*M.^4 + c3*M.^3 + c2*M.^2 + c1*M + c0; 
% Phi_new = alpha.*(1-exp(beta.*(Psi./Psi_th-1)));

a = a3*M.^3 + a2*M.^2 + a1*M + a0;
b = b2*M.^2 + b1*M + b0;
c = c2*M.^2 + c1*M + c0;
Phi_new = a.*(1-exp(Psi.^b+c));

Q = Phi_new*(pi/4)*d_c^2.*U_c;

