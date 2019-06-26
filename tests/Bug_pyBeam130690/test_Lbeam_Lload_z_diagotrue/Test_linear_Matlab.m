clear; clc; close all;
f = 126;
r = 120;

load('ERR_pyBeam_Nastran.mat')

ERR_N_py = (K_N-K_py)./K_N*100;
for i =1:f
    for j = 1:f
        if ( isnan(ERR_N_py(i,j)) || isinf(ERR_N_py(i,j)) )
            ERR_N_py(i,j) = 0;
        end
    end
end

ERR_C_py = (K_C-K_py)./K_C*100;
for i =1:f
    for j = 1:f
        if ( isnan(ERR_C_py(i,j)) || isinf(ERR_C_py(i,j)) )
            ERR_C_py(i,j) = 0;
        end
    end
end

% Full case

R = zeros(f,1);
F = 50000;
DOF = 19*6+3;
R(DOF) = F;

dU_py = K_py\R;
dU_N = K_N\R;
dU_C = K_C\R;






























