% Rotational spring with rigid element
clear; clc;
F = 500000000;
L = 0.3;
E = 70000000000;
ni = 0.3;
G = E/(2*(1+ni));
J = 0.003802666666666667; 
A = L^2*F/G/J;

theta =0.404133297110463;

dy = -L*(1 - cos(theta));
dz = L*sin(theta);

fprintf("dy = %f9. dz = %f9. \n",dy,dz)




