function [sNOxCalc, cNOxCalc, TuCalc, TbCalc, vuCalc, vbCalc, FuCalc, FbCalc, pCylCyc, TCylCyc, FCylCyc, phiCombCylCyc] = NOxPostProcessNew(timeSim,phiCyl,pCyl,TCyl,FCyl,vCyl,mfCyl,QCyl,WCyl,nStroke,varargin)

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

phiCombCyl = mod(phiCyl,nStroke*pi);
phiCombCyl(phiCombCyl>(nStroke/2*pi)) = phiCombCyl(phiCombCyl>(nStroke/2*pi)) - nStroke*pi;

idx = find(phiCombCyl>-pi/9 & phiCombCyl < 1.8*pi/3);
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
% figPosition = [1537.8 501.8 766.4 376.8; ...
%     1537.8 41.8 766.4 376.8; ...
%     2305.8 501.8 766.4 376.8; ...
%     2305.8 41.8 766.4 376.8; ...
%     769.8 41.8 766.4 836.8];
 figPosition = [568.2000  451.4000  444.0000  331.2000; ...
     1013.8    41.8    500.0    328.0; ...
     1012.2    453.0    502.4    328.8; ...
    1515.4    453.8    532.8    328.8; ...
    1515.4    41.8    532.8    328.8];                
                
for i = 1:5
    fig(i) = figure;
    fig(i).Position = figPosition(i,:);
    ax(i) = axes('Parent',fig(i));
    hold(ax(i),'on');
end;    

progress_idx = 0;  
TuCalc = cell(noCyc,1);
TbCalc = cell(noCyc,1);
FuCalc = cell(noCyc,1);
FbCalc = cell(noCyc,1);
vuCalc = cell(noCyc,1);
vbCalc = cell(noCyc,1);
pCylCyc = cell(noCyc,1);
TCylCyc = cell(noCyc,1);
FCylCyc = cell(noCyc,1);
phiCombCylCyc = cell(noCyc,1);
for i = 1:noCyc
    idx = idxStartCyc(i):idxEndCyc(i);
    phiCylCyc = phiCyl(idx);
    phiCombCylCyc{i} = phiCombCyl(idx);
    pCylCyc{i} = pCyl(idx);
    TCylCyc{i} = TCyl(idx);
    FCylCyc{i} = FCyl(idx);
    vCylCyc = vCyl(idx);
    mfCylCyc = mfCyl(idx);
    timeSimCyc = timeSim(idx);
    QCylCyc = QCyl(idx);
    WCylCyc = WCyl(idx);

    [RCylCyc] = GetThdynCombGasZachV1(pCylCyc{i},TCylCyc{i},FCylCyc{i},fs);
    mCyl = pCylCyc{i}.*vCylCyc./ (RCylCyc.*TCylCyc{i});
    nMol = mCyl(end) / (8314.4/RCylCyc(end));  %kmol
    [TuCalc{i},TbCalc{i},FuCalc{i},FbCalc{i},vuCalc{i},vbCalc{i}]= Get2ZoneTemp(phiCylCyc,pCylCyc{i},TCylCyc{i},FCylCyc{i},vCylCyc,mfCylCyc,QCylCyc,FComb,fs,hn);
    plot(ax(1),phiCombCylCyc{i},pCylCyc{i});
    plot(ax(2),phiCombCylCyc{i}(1:end-1),[TuCalc{i} TbCalc{i}]);
    plot(ax(3),phiCombCylCyc{i}(1:end-1),[FuCalc{i} FbCalc{i}]);
    plot(ax(4),phiCombCylCyc{i}(1:end-1),[vuCalc{i} vbCalc{i}]);
    nStep = length(timeSimCyc); 
    idx = 1:(nStep-1);
    NOCalc = zeros(nStep,1);
    NOAcc = 0;
    %Calculation of amount of NO produced per cycle in [mol]
  
    dt = timeSimCyc(2) - timeSimCyc(1);
    dNOFunc =  @(t,NO)dNO(t,NO,timeSimCyc(idx),pCylCyc{i}(idx),TbCalc{i}(idx),FbCalc{i}(idx),vbCalc{i}(idx),fs,xAir,FC,TL,TH);
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [NO] = ode45(dNOFunc,[timeSimCyc(1) timeSimCyc(end-1)],1e-7,options);
    plot(ax(5),interp1(timeSimCyc,phiCombCylCyc{i},NO.x),NO.y);
    NOCyc(i) = NO.y(end);
    cNOCyc(i) = (NOCyc(i))/nMol*1e6; %concentration NO in ppm
    if i - progress_idx > 0.1*noCyc
        progress_idx = progress_idx + 0.1*noCyc;
        fprintf('.');
    end;
end;
%plot(ax(1),cNOCyc)
fprintf('.\n');

idxEndCyc = 0;
noCyc = 1;
for i=1:(length(phiCombCyl)-1)
    if phiCombCyl(i+1) - phiCombCyl(i) < 0  % when the piston at BDC
        noCyc = noCyc + 1;                % count the cycle
        idxEndCyc(noCyc) = i;
    end;
end;
nStep = length(timeSim);
sNOxCalc = zeros(nStep,1);
cNOxCalc = zeros(nStep,1);
sNOxCalc(1:idxEndCyc(1)) = 0;
cNOxCalc(1:idxEndCyc(1)) = 0;


for i=1:length(NOCyc)
    idx = idxEndCyc(i)+1:idxEndCyc(i+1);
    wCyc = (WCyl(idx(end))-WCyl(idx(1)))/3.6e6;          % Work per cycle [kWh]
    sNOxCalc(idx) = NOCyc(i)*46*1e3 / wCyc;                           % Specific NOx [g/kWh]
    cNOxCalc(idx) = cNOCyc(i);
end;
sNOxCalc(idx(end)+1:end) = NOCyc(i)*46*1e3 / wCyc;                           % Specific NOx [g/kWh]
cNOxCalc(idx(end)+1:end) = cNOCyc(i);