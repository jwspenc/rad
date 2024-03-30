
clc
clear, close all

%% Area calc
% MEASUREMENTS UNIT IN [in]
L = 12;
H = 12;
W = 1+5/8;  % For both radiator and fin

% TUBE SHAPE
tuberow = 46;
tubeW = W;
tubeL = 1/16; % [in]

% TUBE AREA (APPROX WATER TUBE TO BE RECTANGULAR)
tubeP = 2*(tubeW+tubeL)/39.3700787; % [in] to [m] perimeter per tube
tubeA = tuberow*2*tubeP*H;  % [in^2] total internal area
tubeA = tubeA/1550; % [m^2] total internal area
tube_oneRowA = tubeA/tuberow;

% FIN SHAPE (APPROX FIN TO BE THE SHAPE OF A SINUSOIDAL WAVE)
finrow = tuberow+1;
finW = W;
finGap = 3/16;  % [in]
finPerin = 16;
freq = finPerin/2;  % [Hz]
omega = 2*pi*freq;  % [rad]

% FIN AREA
fun = @(x) sqrt(1+(3*pi/2*cos(omega*x)).^2);
finL_oneRow = integral(fun, 0, H); % curve length of one period [in]
finL = finL_oneRow*finrow;  % overall fin length [in]
finA = finL*finW;   % [in^2]
finA = finA/1550;   % [m^2] (convert in^2 to m^2 divide by 1550)
fin_oneRowA = finL_oneRow*finW/1550; % [m^2]

% AIR FLOW AREA
airA = (L*H-tuberow*H*tubeL)/1550;  % [m^2] (Total area minus tube longitudinal cross sectional area)

% FIN INTERNAL FLOW CELL AREA
finL_oneT = integral(fun, 0, 1/8);
fin_cellP = finL_oneT+1/8;
fin_cellA = fin_cellP*finW/1550;    % [m^2]
radiator_cellA = fin_cellA*L*freq*finrow;

%% trueAirflow calculated from Fan_Only.mat
trueAirflow = 0.242300722222222; % [m^3/s]

%% HEAT TRANSFER COEFFICIENT CALC FROM NU NUMBER(HTC)
% FIN EXTERNAL FLOW HTC
airVel = trueAirflow/airA; % [m/s]
Re_air = airVel*W*0.0254/20.92e-6; % Laminar
Pr_air = 20.92/29.9;
Nu_air = 0.664*Re_air^(1/2)*Pr_air^(1/3); % Eq. 7.30
h_air = Nu_air*30e-3/W/0.0254; % [W/m^2-K] k_air = 30e-3

% FIN INTERNAL FLOW HTC
Nu_air_internal = 3.73;
Dh_air_internal = 2*1/8*3/16/(1/8+3/16)/39.37;
h_air_internal = Nu_air_internal*30e-3/Dh_air_internal;

% TUBE INTERNAL FLOW HTC
V_tube = 4*0.00378541/60/tubeA; % tube flow velocity
Dh = 4*tubeW*tubeL/1550/tubeP;
Re_water = V_tube*Dh/0.0000003414;
Nu_water = 8.23;
h_water = Nu_water*.674/Dh; % [W/m^2-K] k=0.674 from Table A.4

%% OVERALL HTC FROM THERMAL RESISTANCE
baseA = 2*H*W;
k_Al = 237; % [W/m-K] Table A.1
m = sqrt(2*h_air/k_Al/0.0005); % estimate fin thickness to be 0.0005m (0.5mm)
eta_fin = tanh(m*(finL_oneRow+0.0005/2))/(m*(finL_oneRow+0.0005/2));
eta_0 = 1-fin_oneRowA/baseA*(1-eta_fin);
UA = 1/(1/eta_0/h_air/radiator_cellA+1/h_water/tubeA);

%% HEAT TRANSFER CALCULATION
% CONVERT gal/min TO kg/s
mdot_h = 4*0.00378541/60/0.001034; % 4 * m^3/gal / s/min / m^3/kg

% COOLANT (HOT SIDE) INLET AND OUTLET TEMP
Th_i = 1.260311111111111e+02; % (from Fan_Only.mat)
delTh = 23; % Desired temp drop (used for cp_h estimation)
Th_o = Th_i-delTh;
Th_ave = (Th_i+Th_o)/2+273.15;
cp_h = 4220; % J/kg-K (specific heat at avg temp Table A.6)
Ch = mdot_h*cp_h; % W/K

% AIR (COLD SIDE) INLET AND OUTLET TEMP
Tc_i = 28;
Tc_o = 1.261111111111111e+02; % (from Fan_Only.mat)
Tc_ave = (Tc_i+Tc_o)/2+273.15;
cp_c = 1009; % J/kg-K (specific heat at avg temp Table A.4)
mdot_c = trueAirflow*0.9950; % convert [m^3/s] to [kg/s] (density from Table A.4)
Cc = mdot_c*cp_c;

% HEAT TRANSFER OF BOTH SIDES
q_c = Cc*(Tc_o-Tc_i);
q_h = Ch*(Th_i-Th_o);

%% NTU METHOD
% FIND Cmin
if Ch>Cc
    Cmin = Cc;
    Cmax = Ch;
else
    Cmin = Ch;
    Cmax = Cc;
end

Cr = Cmin/Cmax;
NTU = UA/Cc;
epsilon = 1-exp(1/Cr*NTU^0.22*(exp(-Cr*NTU^0.78)-1));   % double(solve(NTU == -log(1+log(1-epsilon)))); % cross flow both fluid unmixed
q_max = Cc*(Th_i-Tc_i);
q = epsilon*q_max;

del_T = q/Cc
