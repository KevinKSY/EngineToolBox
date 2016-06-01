%% cmapImport.m
%
%   File to input the values of the compressor map. Saves the following values to 'cmap.mat'
%       U2      : corrected speed lines [m/s]
%       U2str   : string containing the values in U2 vector
%       etastr  : string containing the efficiency lines of the compressor map
%       c       : a 3-dimensional matrix with the first column being volume 
%       flow, the second column being pressure ratio and each row giving a
%       different sample. The pages give different speed curves. 
%       Q_eta   : matrix with the volume flow Q_eta(U2,eta) for a given speed
%       line and efficiency
%       eta_double : a repeated vector of efficiency lines = [eta(1) eta(1) eta(2) eta(2) ... eta(end) eta(end)]
%       e       : vector of efficiency samples
%       eQ      : vector of volume flow rate samples corresponding to 'e'
%       eP      : vector of pressure ratio samples corresponding to 'e'
%
% Author: Andreas Torp Karlsen (andreas.t.karlsen@gmail.com)
% Last revised: 02.05.2012

%% Initialization
n = 9; % number of speed curves
m = 6; % number of samples per speed curve
l = 5; % number of efficiency curves

%% The speed curves
%NB: all data values must be distinct to work with interpolation methods

% Three dimensional matrix c(i,j,k) to contain the samples from the speed
% lines of the compressor map
%     i = corr. volume flow [m^3/s]
%     j = pressure ratio [-]
%     k = compressor blade tip speed [m/s]
% Pairs of corrected volume flow and pressure are given for each speed line
c = zeros(m,2,n);

U2 = [250 300 350 400 450 475 500 525 550];  % Speed lines in the compressor map given by compressor blade tip speed [m/s]
U2str = {'U_2=250','U_2=300','U_2=350','U_2=400','U_2=450','U_2=475','U_2=500','U_2=525','U_2=550'}; % String with all the speed lines

%U2 = 250 m/s
c(:,:,1) = [3.35, 1.55;
			3.55, 1.54;
			4.70, 1.52;
			5.50, 1.47;
			6.05, 1.40;
            6.60, 1.28];
        
%U2 = 300 m/s			
c(:,:,2) = [4.60, 1.84;
			4.75, 1.83;
			5.95, 1.80;
			7.00, 1.73;
			7.60, 1.62;
			7.95, 1.40];

%U2 = 350 m/s			
c(:,:,3) = [6.10, 2.22;
			6.11, 2.21;
			7.75, 2.18;
			8.65, 2.07;
			9.30, 1.92;
			9.50, 1.58];
						
%U2 = 400 m/s			
c(:,:,4) = [7.80, 2.73;
			7.90, 2.72;
			9.50, 2.68;
			10.65, 2.54;
			11.20, 2.36;
			11.45, 1.84];
					
%U2 = 450 m/s			
c(:,:,5) = [10.10, 3.42;
			10.15, 3.41;
			11.60, 3.34;
			12.80, 3.12;
			13.25, 2.89;
			13.40, 2.38];
						
%U2 = 475 m/s			
c(:,:,6) = [11.30, 3.82;
			11.40, 3.81;
			12.60, 3.74;
			13.70, 3.49;
			14.05, 3.20;
			14.25, 2.54];
						
%U2 = 500 m/s			
c(:,:,7) = [12.45, 4.24;
			12.55, 4.23;
			13.20, 4.19;
			14.25, 3.93;
			14.55, 3.56;
			14.70, 2.64];	
						
%U2 = 525 m/s			
c(:,:,8) = [13.85, 4.71;
			13.95, 4.67;
			14.05, 4.65;
			14.60, 4.37;
            14.90, 3.78;
			15.05, 2.73];	

%U2 = 550 m/s			
c(:,:,9) = [14.45, 4.85;
			14.60, 4.86;
			14.85, 4.84;
			15.10, 4.52;
			15.30, 4.00;
			15.40, 2.82];	


%% The efficiency curves

% Three dimensional matrix Q(i,j,k) to contain the efficiency curves of the
% compressor map
%     i = speed line number (1...n)
%     j = corrected volume flow rate
%     k = efficiency curve number (1...l)
% For each efficiency curve the two possible volume flow rates are given
% for each speed line it crosses
Q = NaN(n, 2, length(eta));

eta = [0.830 0.820 0.800 0.770 0.730]; % Efficiency curves in the map
etastr = {'0.83', '0.82', '0.80', '0.77', '0.73'}; % String with all the efficiency countours

%eta = 0.830
Q(:,:,1) = [NaN NaN;
            NaN NaN;
            NaN NaN;
            9.45 9.60;
            10.85 12.10;
            11.90 13.10;
            NaN NaN;
            NaN NaN;
            NaN NaN;];
        
%eta = 0.820        
Q(:,:,2) = [NaN NaN;
            NaN NaN;
            6.95 8.10;
            8.30 10.30;
            10.15 12.50;
            11.30 13.45;
            12.90 13.90;
            NaN NaN;
            NaN NaN;];
        
%eta = 0.800        
Q(:,:,3) = [4.25 5.25;
            5.00 6.75;
            6.10 8.60;
            NaN 10.70;
            NaN 12.90;
            NaN 13.80;
            NaN 14.35;
            NaN 14.60;
            NaN NaN;];
        
%eta = 0.770        
Q(:,:,4) = [3.85 5.70;
            NaN 7.20;
            NaN 9.00;
            NaN 11.05;
            NaN 13.20;
            NaN 14.00;
            NaN 14.55;
            NaN 14.90;
            14.60 15.10;];
        
%eta = 0.730     
Q(:,:,5) = [NaN 6.05;
            NaN 7.55;
            NaN 9.25;
            NaN 11.20;
            NaN 13.30;
            NaN 14.10;
            NaN 14.65;
            NaN 14.95;
            NaN 15.30;];
        

% Q_eta (i,j) contains the corrected volume flow for every crossing of
% speed line and efficiency curve.
%    i = speed line (1...n)
%    j = efficiency repeated twice (1...2*l)
% Such that each speed row contains [eta(1) eta(1) eta(2) eta(2) ...
% eta(end) eta(end)] (because there are two possible crossings between each efficiency curve and speed line)
Q_eta = NaN(n,2*length(eta));
eta_double = zeros(2*length(eta),1);

for i=1:n
    for j=1:length(eta)
        Q_eta(i,1+2*(j-1):2+2*(j-1)) = Q(i,:,j);
    end
end

for i=1:length(eta)
    eta_double(1+2*(i-1):2+2*(i-1)) = eta(i);
end

%% Efficiency curves revised

% Pairs (eQ1,eP1) contains the volume flow rate and pressure ratio at an
% arbitrarily point on efficiency curve 1. Must collect enough samples for
% each efficiency curve to be reprodused correctly

eQ1 = [9.50 11.0 12.1 13.1 13.0 11.9 10.8 10.0];
eQ2 = [6.10 8.10 10.25 12.5 13.45 13.90 13.75 12.9 11.3 10.15 8.3 6.95];
eQ3 = [4.5  5.25  6.8 8.65 10.7 12.85 13.8 14.35 14.6 14.05 13.65 6.1  5.0  4.25  4.0];
eQ4 = [4.9 5.7 7.20 8.95 11.05 13.15 14.0 14.55 14.9 15.1 14.6 4.35 3.65];
eQ5 = [6.05 7.55 9.25 11.20 13.3 14.15 14.65 14.95 15.3 15.55];

eP1 = [2.68 2.96 3.27 3.67 4.02 3.8 3.4 3.05];
eP2 = [1.83 2.15 2.6 3.19 3.58 4.08 4.28 4.22 3.82 3.43 2.73 2.21];
eP3 = [1.38 1.49 1.75 2.08 2.53 3.1  3.45 3.86  4.4  4.64  4.63  2.22 1.83 1.53 1.41];
eP4 = [1.34 1.45 1.7 2.02 2.43 2.97 3.27 3.57 3.88 4.52 4.86 1.76 1.45];
eP5 = [1.40 1.63 1.95 2.34 2.82 3.08 3.3 3.5 3.68 4.1];

% eQ, eP and e are vectors containing samples volume flow, pressure ratio and
% efficiency respectively
eQ = [eQ1 eQ2 eQ3 eQ4 eQ5];
eP = [eP1 eP2 eP3 eP4 eP5];

% Allocate the correct efficiency (e) to each index of the triplet (eQ,eP,e)
clear('e') % Just to ensure nothing is collected in variable 'e'
e(1:length(eQ1)) = eta(1);
e(length(e)+1:length(e)+length(eQ2)) = eta(2);
e(length(e)+1:length(e)+length(eQ3)) = eta(3);
e(length(e)+1:length(e)+length(eQ4)) = eta(4);
e(length(e)+1:length(e)+length(eQ5)) = eta(5);

%% Saves the values

save('cmap.mat', 'c', 'U2', 'U2str', 'eta_double', 'Q_eta', 'etastr', 'eQ', 'eP', 'e') 