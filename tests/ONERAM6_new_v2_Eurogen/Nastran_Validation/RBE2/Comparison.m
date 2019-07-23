clc; clear all; close all;
N_nl = [4.056152E-02     -3.172400E-02      1.492862E-01;
    4.056152E-02     -2.106979E-02      1.505786E-01;
    4.056152E-02     -4.770561E-02      1.473475E-01;
    4.071960E-02     -3.751788E-02      1.492862E-01;
    4.040343E-02     -2.593012E-02      1.492862E-01];

py_nl =[0.0406171877987  -0.0317928903361  0.149485266914;
    0.0409804466552  -0.0206963381409  0.152396285513;
    0.0400722900247  -0.0484380169117  0.145118662193;
    0.0406112645021  -0.0375110530265  0.148735502106;
    0.0406231124752  -0.0260747279801  0.150235031634];

py_nl2 =[0.0406151888712  -0.0317898988884  0.149475174739;
    0.040978380444  -0.0206940188178  0.152384572119;
    0.0400703904863  -0.0484340410883  0.145110997087;
    0.0406090848028  -0.0375075611018  0.148725543414;
    0.0406212907351  -0.0260722366515  0.15022480607];

ERR_c1 = (N_nl-py_nl2)./N_nl*100;
%%
N_nl_c2 = [4.285070E-02     -3.307158E-02      1.488857E-01;
    4.285070E-02     -2.274658E-02      1.629326E-01;
    4.285070E-02     -4.855936E-02      1.278149E-01;
    4.456890E-02     -3.958267E-02      1.488857E-01;
    4.113250E-02     -2.656049E-02      1.488857E-01];

py_nl_c2 =[0.0493066373068  -0.0370816334248  0.148410617182;
    0.0553211890985  -0.0366880282485  0.194608910685;
    0.0403712391417  -0.0377398379207  0.0792688640348;
    0.0545170002312  -0.0453680422743  0.146133881904;
    0.0440962607886  -0.028795220581  0.150687357682];

ERR_c2 = (N_nl_c2-py_nl_c2)./N_nl_c2*100;
%%
py_nl_c3 =[-0.000557265823228  -0.00139733860448  0.0589241074485;
    -0.000550809235809  -0.00142928580702  0.0573941936886;
    -0.000566950925416  -0.00134937665216  0.0612189578176;
    -0.000743814639755  -0.00278414565511  0.0588799255415;
    -0.00037072121363  -1.05309945864e-05  0.0589682897976];

N_nl_c3 = [-5.633550E-04     -1.421277E-03      5.937344E-02
    -5.633550E-04     -1.499526E-03      5.781578E-02
    -5.633550E-04     -1.303902E-03      6.170997E-02
    -7.538850E-04     -2.823054E-03      5.937344E-02
    -3.728250E-04     -1.950040E-05      5.937344E-02;];

ERR_c3 = (N_nl_c3-py_nl_c3)./N_nl_c3*100
