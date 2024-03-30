
clc
clear, close all

%% IMPORT DATA (READ MAX & MIN COOLANT TEMPERATURE)
tbl0 = readtable("Howard Wu_WAP 23_Moolah_Race_a_0177.csv");
tbl = table2array(tbl0((1:end-2),(3:end)));
tbl(:,(3:6)) = (tbl(:,(3:6))-32)*5/9;   % Convert temperature [F] to [C]

%% PLOT IMPORTED DATA
% figure;
% plot(tbl(:,1), tbl(:,3))
% ylim([20 130])
% xlabel("Time")
% ylabel("Inlet 1 Temp")
% title("Inlet 1 Temp")
% 
% figure;
% plot(tbl(:,1), tbl(:,4))
% ylim([40 105])
% xlabel("Time")
% ylabel("Inlet 2 Temp")
% title("Inlet 2 Temp")
% 
% figure;
% plot(tbl(:,1), tbl(:,5))
% ylim([30 85])
% xlabel("Time")
% ylabel("Outlet Temp")
% title("Outlet Temp")
% 
% figure;
% plot(tbl(:,1), tbl(:,6))
% ylim([10 130])
% xlabel("Time")
% ylabel("Engine Temp")
% title("Engine Temp")

%% AREA CALCULATION
% MEASUREMENTS UNIT IN [in]
L = 9;
H = 11;

% TUBE MEASUREMENTS
tuberow = 34;
tubeL = 1/16; % [in]

% AIR FLOW AREA
airA = (L*H-tuberow*H*tubeL)/1550;  % [m^2] (Total area minus tube longitudinal cross sectional area)

%% FIND TRUE FLOW RATE (FROM INTERSECTING POINT OF RADIATOR PRESSURE DROP CURVE AND FAN CURVE)
% AIR FLOW RATE CONVERSION (FROM [m/s] TO [m^3/h])
faceVel = [2 4 6];  % Provided by PWR
faceQ = faceVel*airA*3600;
pDrop = [60 145 255];   % Provided by PWR

% FAN CURVE
airQ1 = [1020 930 820 680 450 340 230 150 0];
staticP1 = [0 25 50 75 100 125 150 175 200];
airQ2 = [1850 1670 1490 1090 830 610 460 280 70 0];
staticP2 = (0:50:450);
airQ3 = [1810 1630 1430 1190 860 580 310 130 0];
staticP3 = (0:50:400);
airQ4 = [1830 1630 1430 970 690 540 360 180 60 0];
staticP4 = (0:50:450);
airQ5 = [1750 1560 1380 1120 810 540 300 210 0];
staticP5 = (0:50:400);
airQ6 = [2720 2570 2420 2270 2100 1830 1570 1030 840 630 410 270 0];
staticP6 = (0:50:600);

f2 = figure;
plot(airQ1, staticP1)   % current fan
hold on
plot(airQ2, staticP2)   % VA17-AP70-/LL-39A [10in]
plot(airQ3, staticP3)   % VA17-AP70-/LL-39S [10in]
plot(airQ4, staticP4)   % VA17-AP70-/LL-51A [10in]
plot(airQ5, staticP5)   % VA15-AP70/LL-39S [10in] (chosen)
plot(airQ6, staticP6)   % VA03-AP90-/LL-68A [11in]
plot(faceQ, pDrop)
plot([faceQ(end) 1500],[pDrop(end) 380],"-.")
grid on
% xline(495.9581)
% xline(872.2826)
xline(855.0391)
xline(1300)
text(857,580,'855.0391')
text(1300,580,'1300')
hold off
xlabel('Airflow [m^3/h]')
ylabel('Static pressure [Pa]')
legend('Current Fan', 'VA17-AP70-/LL-39A', 'VA17-AP70-/LL-39S',...
    'VA17-AP70-/LL-51A','VA15-AP70/LL-39S','VA03-AP90-/LL-68A',...
    'Rad pressure loss','Extrapolated loss')

% exportgraphics(f2,"Rad and Fan curve.png",'Resolution',600)

% FIND INTERSECTING POINT IN PLOT AS TRUE AIR FLOW RATE
trueAirflow = 855.0391/3600; % [m^3/s] (divide by 3600 [m^3/h] to [m^3/s]) 495.9581 872.2826 855.0391

%% HEAT TRANSFER CALCULATION (Eq. 11.6b & 11.7b)
% AIR (COLD SIDE) INLET AND OUTLET TEMP
Tc_i = 28;
Tc_o = max(tbl(:,6));
Tc_ave = (Tc_i+Tc_o)/2+273.15;
cp_c = 1009; % [J/kg-K] (specific heat at avg temp Table A.4)
mdot_c = trueAirflow*0.9950; % [kg/s] Convert Q[m^3/s] to mdot[kg/s] (density from Table A.4)
Cc = mdot_c*cp_c;

% COOLANT (HOT SIDE) INLET AND OUTLET TEMP
Th_i = max(tbl(:,3));
cp_h = 4220; % [J/kg-K] Best guess (specific heat from Table A.6)
mdot_h = 4*0.00378541/60/0.001034; % [kg/s] Water pump flow rate: 4gpm * [m^3/gal] / [s/min] / [m^3/kg]
Ch = mdot_h*cp_h;

% HEAT TRANSFER OF COLD SIDE
q_c = Cc*(Tc_o-Tc_i);   % [W]

% FIND POSSIBLE TEMP DROP OF HOT SIDE
delTh = q_c/Ch
Th_o = Th_i-delTh;  % Radiator outlet temperature

% CAN USE THE FOLLOWING CODE TO CHECK WATER PROPERTY
% Th_ave = (Th_i+Th_o)/2+273.15;