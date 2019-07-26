clear; clc; close all;
% 
% % Finite difference evaluation   % Nominal Load 100000
% % DeltaF
% 
V_nominal =21.51333257972762424970;


dv = 21.51342713648649862535; dP = 1.0e-0;

dv_dP = (dv - V_nominal)/dP                

dv_dP_AD = 9.455735874592986e-05;

ERR = (dv_dP - dv_dP_AD)/dv_dP*100

%% Young modulus

dv = 21.51333257959254652292; dP = 1.0e-0;

dv_dP = (dv - V_nominal)/dP                

dv_dP_AD = -1.3508194106561384e-10;

ERR = (dv_dP - dv_dP_AD)/dv_dP*100








