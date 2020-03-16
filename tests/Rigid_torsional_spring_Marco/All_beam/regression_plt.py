#!/usr/bin/env python3
#
# pyBeam, a Beam Solver
#
# Copyright (C) 2018 Ruben Sanchez, Rocco Bombardieri, Rauno Cavallaro
# 
# Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
#             Rocco Bombardieri, Rauno Cavallaro (Carlos III University Madrid)
#
# This file is part of pyBeam.
#
# pyBeam is free software: you can redistribute it and/or
# modify it under the terms of the GNU Affero General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# pyBeam is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero
# General Public License along with pyBeam.
# If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys, os
from pyBeamLib import pyBeamSolver

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))

beam = pyBeamSolver(file_dir ,'config_NL.cfg')

beam.SetLoads(2 , 0, 0 , 500000000) #
#beam.SetLoads(19 , 0, 0 , 7.8e-02)
#beam.SetLoads(19 , 0.78005, 0 , 0)

beam.Run()

beam.coordinate_Y1 = beam.coordinate_Y; beam.coordinate_X1 = beam.coordinate_X; beam.coordinate_Z1 = beam.coordinate_Z


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create cubic bounding box to simulate equal aspect ratio
max_range = np.array([np.amax(beam.coordinate_X0) - np.amin(beam.coordinate_X0), np.amax(beam.coordinate_Y0) - np.amin(beam.coordinate_Y0),
                      np.amax(beam.coordinate_Z0) - np.amin(beam.coordinate_Z0)]).max()
Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (
            np.amax(beam.coordinate_X0) + np.amin(beam.coordinate_X0))
Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (
            np.amax(beam.coordinate_Y0) + np.amin(beam.coordinate_Y0))
Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (
            np.amax(beam.coordinate_Z0) + np.amin(beam.coordinate_Z0))
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
    plt.plot([xb], [yb], [zb], 'w')

#plt.plot(beam.coordinate_X2[0:20], beam.coordinate_Y2[0:20], beam.coordinate_Z2[0:20],'r')
plt.plot(beam.coordinate_X1[0:20], beam.coordinate_Y1[0:20], beam.coordinate_Z1[0:20],'b')
plt.plot(beam.coordinate_X0[0:20], beam.coordinate_Y0[0:20], beam.coordinate_Z0[0:20],'g')
rigid = 1

'''
for i in range(1,1+rigid):
    node_i = i
    node_j = i+1
    #plt.plot([ beam.coordinate_X2[node_i],beam.coordinate_X2[node_j] ], [ beam.coordinate_Y2[node_i],beam.coordinate_Y2[node_j] ], [ beam.coordinate_Z2[node_i],beam.coordinate_Z2[node_j] ],'r')
    plt.plot([ beam.coordinate_X1[node_i],beam.coordinate_X1[node_j] ], [ beam.coordinate_Y1[node_i],beam.coordinate_Y1[node_j] ], [ beam.coordinate_Z1[node_i],beam.coordinate_Z1[node_j] ],'b')
    plt.plot([ beam.coordinate_X0[node_i],beam.coordinate_X0[node_j] ], [ beam.coordinate_Y0[node_i],beam.coordinate_Y0[node_j] ], [ beam.coordinate_Z0[node_i],beam.coordinate_Z0[node_j] ],'g')
'''


plt.xlabel('X')
plt.ylabel('Y')
plt.show()

for i in range(beam.nPoint):
    print('{} {}  {}  {} '.format(i, beam.displacement_X[i],beam.displacement_Y[i],beam.displacement_Z[i]) )
