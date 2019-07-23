clear; clc; close all;
% 
% % Finite difference evaluation   % Nominal Load 100000  % Nominal YOUNG MODULUS 70000000000
% % DeltaE
% 
V_nominal = 21.51333254587530419144;

dv = 21.51333254587530419144; dP = 1.0e-8;

dv_dP_1 = (dv - V_nominal)/dP

dv = 21.51333254587530419144; dP = 1.0e-6;

dv_dP_2 = (dv - V_nominal)/dP


dv = 21.51333248026568867317; dP = 1.0e-4;

dv_dP_3 = (dv - V_nominal)/dP


dv = 21.51333247360780376312; dP = 1.0e-3;

dv_dP_4 = (dv - V_nominal)/dP

dv = 21.51333249227407762305; dP = 1.0e-2;

dv_dP_5 = (dv - V_nominal)/dP

dv = 21.51333246149479805354; dP = 1.0e-1;

dv_dP_6 = (dv - V_nominal)/dP

dv = 21.51333252495046721720; dP = 1.0e-0;

dv_dP_7 = (dv - V_nominal)/dP

dv = 21.51333248186155344683; dP = 1.0e+2;

dv_dP_8 = (dv - V_nominal)/dP

dv = 21.51333236687237970841; dP = 1.0e+3;

dv_dP_9 = (dv - V_nominal)/dP

dv = 21.51333112589327356545; dP = 1.0e+4;

dv_dP_10 = (dv - V_nominal)/dP

dv = 21.51331911355023152055; dP = 1.0e+5;

dv_dP_11 = (dv - V_nominal)/dP

dv = 21.51319741933159335190; dP = 1.0e+6;

dv_dP_12 = (dv - V_nominal)/dP

dv = 21.51198175789615163467; dP = 1.0e+7;

dv_dP_13 = (dv - V_nominal)/dP    %% Winner


dv_dP_AD = -1.3508191358565246e-10;

ERR = ( dv_dP_13 - dv_dP_AD)/dv_dP_13*100















