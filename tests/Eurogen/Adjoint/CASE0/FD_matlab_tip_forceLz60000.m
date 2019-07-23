clear; clc; close all;
% 
% % Finite difference evaluation   % Nominal Load 60000
% % DeltaF
% 
V_nominal = 16.45455992503018904927;

dv = 16.45457598153281253417; dP = 1.0e-1;

dv_dP_1 = (dv - V_nominal)/dP

dv = 16.45472934010101084823; dP = 1.0e-0;

dv_dP_2 = (dv - V_nominal)/dP

dv = 16.45624600694933192813; dP = 1.0e+1;

dv_dP_3 = (dv - V_nominal)/dP

dv = 16.47138454191162182383; dP = 1.0e+2;

dv_dP_4 = (dv - V_nominal)/dP

dv = 16.62161475318926306954; dP = 1.0e+3;

dv_dP_5 = (dv - V_nominal)/dP

dv = 16.46297531139440550874; dP = 5.0e+1;    

dv_dP_6 = (dv - V_nominal)/dP                 %% Winner

dv_dP_AD = 0.00016840030716481512;

ERR = (dv_dP_6 - dv_dP_AD)/dv_dP_6*100

% %% Linear case
% 
% L = 30;
% 
% P = 100000; E = 70000000000; I = 0.0003533333333333333;
% 
% v= P*L^3/3/E/I
% 
% dv_dP_analytic = L^3/3/E/I
% 
% V_nominal = 36.38814016182118393772;
% 
% %% FINITE DIFFERENCES
% dv = 36.38814016182122657028; dP = 1.0e-10;
% 
% dv_dP_1 = (dv - V_nominal)/dP
% 
% 
% dv = 36.38814016182156052537; dP = 1.0e-9;
% 
% dv_dP_2 = (dv - V_nominal)/dP
% 
% dv = 36.38814016182483612738; dP = 1.0e-8;
% 
% dv_dP_3 = (dv - V_nominal)/dP
% 
% dv = 36.38814016185756372579; dP = 1.0e-7;
% 
% dv_dP_4 = (dv - V_nominal)/dP
% 
% dv = 36.38814016218507418898; dP = 1.0e-6;
% 
% dv_dP_5 = (dv - V_nominal)/dP
% 
% dv = 36.38814016546000829067; dP = 1.0e-5;
% 
% dv_dP_5 = (dv - V_nominal)/dP





