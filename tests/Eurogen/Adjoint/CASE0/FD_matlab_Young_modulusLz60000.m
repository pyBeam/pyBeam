clear; clc; close all;
% 
% % Finite difference evaluation   % Nominal Load 60000  % Nominal YOUNG MODULUS 70000000000
% % DeltaE
% 
V_nominal = 16.45455992503018904927;

dv = 16.45441753963913456005; dP = 1.0e+6;

dv_dP_1 = (dv - V_nominal)/dP

dv = 16.45311884117913336922; dP = 1.0e+7;

dv_dP_2 = (dv - V_nominal)/dP

dv = 16.44013645603416762242; dP = 1.0e+8;

dv_dP_3 = (dv - V_nominal)/dP


dv_dP_AD = -1.44343120426939e-10;

ERR = ( dv_dP_3 - dv_dP_AD)/dv_dP_3*100















