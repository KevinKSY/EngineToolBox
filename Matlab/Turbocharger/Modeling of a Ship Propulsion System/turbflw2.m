function Q = turbflw2(k,PR,k_e)
% Q = TURBFLW2(k,PR)
% Calculates turbine flow [m^3/(s sqrtK)]
%   Q = k*Psi(PR) or Q = (k(1)/PR + k(2))*Psi
%
%   Data:
%       k   : parameter(s) to be fitted
%       PR  : pressure expansion ratio [-]
%       k_e : ratio of specific heats for exhaust
%   Output:
%       Q   : turbine flow [m^3/(s sqrtK)]
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 05/2012

Psi = zeros(1,length(PR));
for i=1:length(Psi)
    if PR(i) > ((k_e+1)/2)^(k_e/(k_e-1))
        Psi(i) = sqrt(k_e*(2/(k_e+1))^((k_e+1)/(k_e-1)));
    else
        Psi(i) = (1/PR(i))^(1/k_e)*sqrt( (2*k_e)/(k_e-1)*( 1 - (1/PR(i))^((k_e-1)/k_e) ) );
    end
end
    
%Q = k(1).*Psi;
Q = (k(1)./PR + k(2)).*Psi;

end