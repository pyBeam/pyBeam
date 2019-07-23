clear; clc; close all;
% 
% % Finite difference evaluation   % Nominal Load 100000
% % DeltaF
% 
V_nominal = 21.51333254587530419144;

dv = 21.51333247228767575621; dP = 1.0e-8;

dv_dP_1 = (dv - V_nominal)/dP

dv = 21.51333248228514349876; dP = 1.0e-7;

dv_dP_2 = (dv - V_nominal)/dP

dv = 21.51333258007733917339; dP = 1.0e-6;

dv_dP_3 = (dv - V_nominal)/dP

dv = 21.51333253318999894077; dP = 1.0e-5;

dv_dP_4 = (dv - V_nominal)/dP

dv = 21.51333250149011178110; dP = 1.0e-4;

dv_dP_5 = (dv - V_nominal)/dP

dv = 21.51333259712772516536; dP = 1.0e-3;

dv_dP_6 = (dv - V_nominal)/dP

dv = 21.51333342191078301653; dP = 1.0e-2;

dv_dP_7 = (dv - V_nominal)/dP

dv = 21.51334193818130202658; dP = 1.0e-1;

dv_dP_8 = (dv - V_nominal)/dP

dv = 21.51342709474404912839; dP = 1.0e-0;

dv_dP_9 = (dv - V_nominal)/dP                %%% Winner

dv = 21.51427795407120413529; dP = 1.0e+1;

dv_dP_10 = (dv - V_nominal)/dP

dv = 21.52278229195096415083; dP = 1.0e+2;

dv_dP_11 = (dv - V_nominal)/dP

dv_dP_AD = 9.455733950996038e-05;

ERR = (dv_dP_9 - dv_dP_AD)/dv_dP_9*100


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





