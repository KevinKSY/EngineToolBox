function Q = turbflw(k,PR)
% Q = TURBFLW(k,PR)
% Calculates turbine flow [m^3/(s sqrtK)]
%   Q = k(1).*sqrt(1-PR.^k(2)); 
%
%   Data:
%       k   : vector of parameters to be fitted (1x2)
%       PR  : pressure expansion ratio [-]
%   Output:
%       Q   : turbine flow [m^3/(s sqrtK)]
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Revised: 05/2012

Q = k(1).*sqrt(1-PR.^k(2));

end