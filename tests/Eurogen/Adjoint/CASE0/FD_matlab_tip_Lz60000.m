clear; clc; close all;
% 
% % Finite difference evaluation   % Nominal Load 60000
% % DeltaF
% 
V_nominal = 16.45454581459997456250;

%% Force

dv = 16.462962234119175; dP = 5.0e+1;    

dv_dP_6 = (dv - V_nominal)/dP                

dv_dP_AD = 0.000168396253822679;

ERR = (dv_dP_6 - dv_dP_AD)/dv_dP_6*100


%% Young Modulus

dv = 16.45454581445567754372; dP = 1.0e+0;

dv_dP_3 = (dv - V_nominal)/dP


dv_dP_AD = -1.4433964613372497e-10;

ERR_E = ( dv_dP_3 - dv_dP_AD)/dv_dP_3*100

