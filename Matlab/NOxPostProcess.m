function [NOCalc, sNOxCalc, sumNOx] = NOxPostProcess(timeSim,phiCyl,pCyl,TCyl,FCyl,vCyl,mfCyl,QCyl,WCyl,varargin)

if isempty(varargin)
    FComb = 1;                  %Fuel air equivalent ratio of combustion
    fs = 0.0683;                %fs
    xAir = [0;0;0;0;0;0;0;0.2095;0;0.0003;0.7809;0.0093]; %Air composition in mole fraction
    FC = [1;1.8;0;0];           %Fuel composition as in [m;n;p;q] where CmHnOpNq
    TL = 1500;                  %Tempreature low limit for transition from complete combustion to dissociation
    TH = 1700;                  %Temperature high limit for transition from complete combustion to dissociation
    hn = 42707000;
else
    if isempty(varargin{1})
        FComb = 1;
    else 
        FComb = varargin{1};
    end;
    if isempty(varargin{2})
        fs = 0.0683;
    else
        fs = varargin{2};
    end;
    if isempty(varargin{3})
        xAir = [0;0;0;0;0;0;0;0.2095;0;0.0003;0.7809;0.0093];
    else
        xAir = varargin{3};
    end;    
    if isempty(varargin{4})
        FC = [1;1.8;0;0];
    else
        FC = varargin{4};
    end;
    if isempty(varargin{5})
        TL = 1500; TH = 1700;
    else
        TL = varargin{5}(1);
        TH = varargin{5}(2);
    end;    
    if isempty(varargin{6})
        hn = 42707000;
    else
        hn = varargin{6}*1e6;
    end;        
end;


[RCyl] = GetThdynCombGasZach(pCyl,TCyl,FCyl);
mCyl = pCyl.*vCyl./ (RCyl.*TCyl);


[TuCalc,TbCalc,FuCalc,FbCalc,vuCalc,vbCalc]= Get2ZoneTemp(phiCyl,pCyl,TCyl,FCyl,vCyl,mfCyl,QCyl,FComb,fs,hn);

timeSim = timeSim(1:end-1);
pCyl = pCyl(1:end-1);
mCyl = mCyl(1:end-1);
nStep = length(timeSim);
NOCalc = zeros(nStep,1);
sumNOx = zeros(nStep,1);

phiCombCyl = phiCyl - 2*pi*(floor(phiCyl/(2*pi)));
phiCombCyl(phiCombCyl>pi) = phiCombCyl(phiCombCyl>pi) - 2*pi;

idx = find(phiCombCyl>-pi/9 & phiCombCyl < pi/2);
noIdx = length(idx);
noCyc = 1;
idxStartCyc(noCyc) = idx(1);
for i=1:noIdx-1
    if idx(i+1) - idx(i) ~= 1;
        idxEndCyc(noCyc) = idx(i);
        noCyc = noCyc + 1;
        idxStartCyc(noCyc) = idx(i+1);
    end;
end;
if noCyc ~= length(idxEndCyc) 
    noCyc = noCyc - 1;
end;


NOAcc = 0;
%Calculation of amount of NO produced per cycle in [mol]
progress_idx = 0;
for i = 1:noCyc
    idx = idxStartCyc(i):idxEndCyc(i);
    dNOFunc =  @(t,NO)dNO(t,NO,timeSim(idx),mCyl(idx),pCyl(idx),TbCalc(idx),FbCalc(idx),xAir,FC,TL,TH);
    ts = timeSim(idx(1));
    tf = timeSim(idx(end));
    [t,NO] = ode45(dNOFunc,[ts tf],0);
    NOAcc = NOAcc + max(NO)*30;                 % NOx produced 
    sumNOx(idx) = NOAcc;
    NOCalc(idx) = interp1(t,NO,timeSim(idx));
    if i - progress_idx > 0.1*noCyc
        progress_idx = progress_idx + 0.1*noCyc;
        fprintf('.');
    end;
end;
fprintf('.\n');

idxEndCyc(1) = 0;
noCyc = 1;
for i=1:nStep-1
    if phiCombCyl(i+1) - phiCombCyl(i) < 0  % when the piston at BDC
        noCyc = noCyc + 1;                % count the cycle
        idxEndCyc(noCyc) = i;
    end;
end;
sNOxCalc = zeros(nStep,1);
for i=1:noCyc-1
    idx = idxEndCyc(i)+1:idxEndCyc(i+1);
    NO = max(NOCalc(idx))*30;                   % NOx produced per cycle [g]
    wCyc = WCyl(idx(end))-WCyl(idx(1));          % Work per cycle [kWh]
    sNOx = NO / wCyc;                           % Specific NOx [g/kWh]
    sNOxCalc(idx) = sNOx;
end;