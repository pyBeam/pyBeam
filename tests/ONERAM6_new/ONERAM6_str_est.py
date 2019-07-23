#!/usr/bin/env python

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                              %
# BEAM configuration file                                                      %
# Case description: ONERAM6 routine for structural parameer estimation         %
# Author: _Rocco Bombardieri ________________________________________________  %
# Institution: _UC3M_________________________________________________________  %
# Date: _11/06/19_                                                             %
# File Version                                                                 %
#                                                                              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import sys, os
import numpy as np
from math import *

# Set 3 stations along the span
# We can consider the wingbox section length to be 3/4 of the chord and the height to be 4/5 of the thickness

# Panel thickness is assigned so that for a given vertical load applied to the wing tip the vertical deflection is 20% of the span
t = 0.0005   #[m]


# y = 0
c1 = 0.806   #[m]
t1 = 0.078   #[m]
b1 = 3*c1/4
h1 = 4*t1/5

# y = 0.6 [m]
c2 = 0.657   #[m]
t2 = 0.064   #[m]
#
b2 = 3*c2/4
h2 = 4*t2/5

# y = 1.196 m
c3 = 0.454   #[m]
t3 = 0.044   #[m]
#
b3 = 3*c3/4
h3 = 4*t3/5


# Evaluation of the section properties
A1 = 2*(b1+h1)*t
A2 = 2*(b2+h2)*t
A3 = 2*(b3+h3)*t

print("b1 = {}, h1 = {}".format(b1,h1))
print("b2 = {}, h2 = {}".format(b2,h2))
print("b3 = {}, h3 = {}".format(b3,h3))

print("A1 = {}, A2 = {}, A3 = {}".format(A1,A2,A3))

Ix1 = 2*h1**3*t/12 + 2*t**3*b1/12 + 2*b1*t*(h1/2)**2
Ix2 = 2*h2**3*t/12 + 2*t**3*b2/12 + 2*b2*t*(h2/2)**2
Ix3 = 2*h3**3*t/12 + 2*t**3*b3/12 + 2*b3*t*(h3/2)**2

print("Ix1 = {}, Ix2 = {}, Ix3 = {}".format(Ix1,Ix2,Ix3))

Iz1 = 2*b1**3*t/12 + 2*t**3*h1/12 + 2*h1*t*(b1/2)**2
Iz2 = 2*b2**3*t/12 + 2*t**3*h2/12 + 2*h2*t*(b2/2)**2
Iz3 = 2*b3**3*t/12 + 2*t**3*h3/12 + 2*h3*t*(b3/2)**2

print("Iz1 = {}, Iz2 = {}, Iz3 = {}".format(Iz1,Iz2,Iz3))

Jt1 = 2*t**2*(b1-t)**2*(h1-t)**2/(h1*t+b1*t-2*t**2)
Jt2 = 2*t**2*(b2-t)**2*(h2-t)**2/(h2*t+b2*t-2*t**2)
Jt3 = 2*t**2*(b3-t)**2*(h3-t)**2/(h3*t+b3*t-2*t**2)



print("Jt1 = {}, Jt2 = {}, Jy3 = {}".format(Jt1,Jt2,Jt3))







