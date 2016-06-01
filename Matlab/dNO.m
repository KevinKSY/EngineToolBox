function dNO = dNOBowman(t,NO,timeSim,pCyl,TCyl,FCyl,vCyl,fs,varargin)
% The function calculates the derivative of NO formation by thermal NOx 
% using Zedovich's extended mechanism.
% Input (Variable)
%   t : simulation time
%   NO : NO formed          [mol]
% Input (Parameters)
%   tArray : time array for given cylinder data  [s]
%   pCyl : Time-series of pressure of the volume [Pa]
%   TCyl : Time-series of temperature of the volume [K]
%   FCyl : Time-series of fuel-air equivalent ratio of the volume 
% Input (Optional)
%   xAir : Air composition in mole fraction 
%           default([0;0;0;0;0;0;0;0.2095;0;0.0003;0.7809;0.0093])
%   FC : Fuel composition as in [m;n;p;q] where CmHnOpNq
%           default([1;1.8;0;0])          
%   TL : Tempreature low limit for transition from complete combustion to
%        dissociation (default 1500)                  
%   TH : Tempreature high limit for transition from complete combustion to
%        dissociation (default 1700)
% Output
%   dNO : Time derivative of the NO formation [mol/s]

if isempty(varargin)
    xAir = [0;0;0;0;0;0;0;0.2095;0;0.0003;0.7809;0.0093];
                            %Composition of the air (mole fraction) 
                            % [xH; xO; xN; xH2; xOH; xCO; xNO; xO2; xH2O;
                            %     xCO2; xN2; xAr];
    FC = [1;1.8;0;0];        %Composition of the fuel (CmHnOoNp)
    T_L = 1500;
    T_H = 1700;
else
    if isempty(varargin{1})
        xAir = [0;0;0;0;0;0;0;0.2095;0;0.0003;0.7809;0.0093];
    else
        xAir = varargin{1};        
    end;
    if isempty(varargin{2})
        FC = [1;1.8;0;0];        %Composition of the fuel (CmHnOoNp)
    else
        FC = varargin{2};
    end;
    if isempty(varargin{3})
        T_L = 1500;
    else
        T_L = varargin{3};
    end;
    if isempty(varargin{4})
        T_H = 1700;
    else
        T_H = varargin{4};        
    end;
end;
idx = 1;
while t > timeSim(idx+1)
    idx = idx + 1;
end;
if idx >= length(pCyl)
    idx = idx - 1;
end;
didx = (t - timeSim(idx))/(timeSim(idx+1)-timeSim(idx));
p = pCyl(idx) + didx*(pCyl(idx+1) - pCyl(idx));
T = TCyl(idx) + didx*(TCyl(idx+1) - TCyl(idx));
F = FCyl(idx) + didx*(FCyl(idx+1) - FCyl(idx));
v = (vCyl(idx) + didx*(vCyl(idx+1) - vCyl(idx)));                          %m3
[R,~,~,~] = GetThdynCombGasZachV1(p,T,F,fs);
MW = 8314.4621/R;                      %kg/kmol
m = p*v/(R*T);                          %kg
noMol = m/MW;                               %kmol

%v = m*R*T/p*1e6;                            %volume in cm3
x = GetCompCombGas(p,T,F,T_L,T_H,FC,xAir);  %Composition in mole fraction
N = x*noMol;                                %Composition in kmol
c = N/v;                                  % concentration kmol/m3
k1 = 7.6e10*exp(-38000/T);
k2 = 6.4e6*T*exp(-3150/T);
k3 = 4.1e10;
R1 = k1*c(2)*c(11);
R2 = k2*c(3)*c(8);
R3 = k3*c(3)*c(5);
cNO = NO / v;
if R2+R3 == 0
    dcNO = 0;
else
    dcNO = 2*R1*(1-(cNO/c(7))^2) / (1+cNO/c(7)*R1/(R2+R3));
end;
dNO = dcNO*v;



