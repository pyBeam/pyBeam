#!/usr/bin/env python
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

iNode = 59

beam = pyBeamSolver(file_dir ,'config_NL.cfg')

beam.SetLoads(        0 ,  0.00013,  0.00004 ,  0.00114 )
beam.SetLoads(        1 ,  0.00010,  0.00003 ,  0.00122 )
beam.SetLoads(        2 ,  0.00009,  0.00002 ,  0.00129 )
beam.SetLoads(        3 ,  0.00007,  0.00001 ,  0.00135 )
beam.SetLoads(        4 ,  0.00005,  0.00000 ,  0.00139 )
beam.SetLoads(        5 ,  0.00004,  0.00000 ,  0.00142 )
beam.SetLoads(        6 ,  0.00002,  0.00000 ,  0.00143 )
beam.SetLoads(        7 ,  0.00001,  0.00000 ,  0.00144 )
beam.SetLoads(        8 , -0.00000,  0.00001 ,  0.00143 )
beam.SetLoads(        9 , -0.00001,  0.00002 ,  0.00141 )
beam.SetLoads(       10 , -0.00002,  0.00003 ,  0.00137 )
beam.SetLoads(       11 , -0.00003,  0.00004 ,  0.00132 )
beam.SetLoads(       12 , -0.00004,  0.00006 ,  0.00126 )
beam.SetLoads(       13 , -0.00004,  0.00008 ,  0.00119 )
beam.SetLoads(       14 , -0.00005,  0.00010 ,  0.00110 )
beam.SetLoads(       15 , -0.00005,  0.00012 ,  0.00100 )
beam.SetLoads(       16 , -0.00005,  0.00015 ,  0.00089 )
beam.SetLoads(       17 , -0.00005,  0.00018 ,  0.00076 )
beam.SetLoads(       18 , -0.00005,  0.00021 ,  0.00063 )
beam.SetLoads(       19 , -0.00005,  0.00024 ,  0.00048 )
beam.SetLoads(       20 ,  0.00014, -0.00006 ,  0.00153 )
beam.SetLoads(       21 ,  0.00011, -0.00006 ,  0.00167 )
beam.SetLoads(       22 ,  0.00009, -0.00006 ,  0.00179 )
beam.SetLoads(       23 ,  0.00006, -0.00005 ,  0.00189 )
beam.SetLoads(       24 ,  0.00004, -0.00005 ,  0.00198 )
beam.SetLoads(       25 ,  0.00002, -0.00004 ,  0.00206 )
beam.SetLoads(       26 , -0.00000, -0.00003 ,  0.00212 )
beam.SetLoads(       27 , -0.00002, -0.00002 ,  0.00216 )
beam.SetLoads(       28 , -0.00004, -0.00000 ,  0.00218 )
beam.SetLoads(       29 , -0.00005,  0.00001 ,  0.00219 )
beam.SetLoads(       30 , -0.00007,  0.00003 ,  0.00219 )
beam.SetLoads(       31 , -0.00008,  0.00005 ,  0.00216 )
beam.SetLoads(       32 , -0.00009,  0.00008 ,  0.00213 )
beam.SetLoads(       33 , -0.00010,  0.00010 ,  0.00207 )
beam.SetLoads(       34 , -0.00011,  0.00013 ,  0.00200 )
beam.SetLoads(       35 , -0.00011,  0.00016 ,  0.00192 )
beam.SetLoads(       36 , -0.00012,  0.00019 ,  0.00181 )
beam.SetLoads(       37 , -0.00012,  0.00022 ,  0.00170 )
beam.SetLoads(       38 , -0.00012,  0.00026 ,  0.00156 )
beam.SetLoads(       39 , -0.00012,  0.00030 ,  0.00141 )
beam.SetLoads(       40 , -0.00002,  0.00008 ,  0.00032 )
beam.SetLoads(       41 , -0.00002,  0.00005 ,  0.00033 )
beam.SetLoads(       42 , -0.00002,  0.00003 ,  0.00032 )
beam.SetLoads(       43 , -0.00003,  0.00001 ,  0.00031 )
beam.SetLoads(       44 , -0.00003, -0.00001 ,  0.00030 )
beam.SetLoads(       45 , -0.00003, -0.00002 ,  0.00027 )
beam.SetLoads(       46 , -0.00003, -0.00003 ,  0.00023 )
beam.SetLoads(       47 , -0.00003, -0.00004 ,  0.00019 )
beam.SetLoads(       48 , -0.00003, -0.00005 ,  0.00014 )
beam.SetLoads(       49 , -0.00003, -0.00005 ,  0.00007 )
beam.SetLoads(       50 , -0.00003, -0.00004 ,  0.00001 )
beam.SetLoads(       51 , -0.00002, -0.00004 , -0.00007 )
beam.SetLoads(       52 , -0.00002, -0.00003 , -0.00016 )
beam.SetLoads(       53 , -0.00002, -0.00002 , -0.00025 )
beam.SetLoads(       54 , -0.00001, -0.00000 , -0.00036 )
beam.SetLoads(       55 , -0.00001,  0.00002 , -0.00047 )
beam.SetLoads(       56 , -0.00000,  0.00004 , -0.00059 )
beam.SetLoads(       57 ,  0.00000,  0.00007 , -0.00072 )
beam.SetLoads(       58 ,  0.00001,  0.00009 , -0.00086 )
beam.SetLoads(       59 ,  0.00002,  0.00013 , -0.00100 )
beam.SetLoads(       60 , -0.00000,  0.00015 ,  0.00628 )
beam.SetLoads(       61 , -0.00003,  0.00014 ,  0.00627 )
beam.SetLoads(       62 , -0.00004,  0.00013 ,  0.00624 )
beam.SetLoads(       63 , -0.00006,  0.00012 ,  0.00621 )
beam.SetLoads(       64 , -0.00008,  0.00011 ,  0.00615 )
beam.SetLoads(       65 , -0.00009,  0.00011 ,  0.00609 )
beam.SetLoads(       66 , -0.00011,  0.00011 ,  0.00601 )
beam.SetLoads(       67 , -0.00012,  0.00011 ,  0.00591 )
beam.SetLoads(       68 , -0.00013,  0.00012 ,  0.00580 )
beam.SetLoads(       69 , -0.00014,  0.00012 ,  0.00568 )
beam.SetLoads(       70 , -0.00015,  0.00013 ,  0.00554 )
beam.SetLoads(       71 , -0.00016,  0.00014 ,  0.00539 )
beam.SetLoads(       72 , -0.00016,  0.00016 ,  0.00522 )
beam.SetLoads(       73 , -0.00016,  0.00017 ,  0.00504 )
beam.SetLoads(       74 , -0.00017,  0.00019 ,  0.00485 )
beam.SetLoads(       75 , -0.00017,  0.00022 ,  0.00464 )
beam.SetLoads(       76 , -0.00017,  0.00024 ,  0.00442 )
beam.SetLoads(       77 , -0.00016,  0.00027 ,  0.00418 )
beam.SetLoads(       78 , -0.00016,  0.00030 ,  0.00393 )
beam.SetLoads(       79 , -0.00016,  0.00033 ,  0.00367 )
beam.SetLoads(       80 ,  0.00026, -0.00007 , -0.00399 )
beam.SetLoads(       81 ,  0.00024, -0.00008 , -0.00382 )
beam.SetLoads(       82 ,  0.00022, -0.00009 , -0.00366 )
beam.SetLoads(       83 ,  0.00020, -0.00010 , -0.00351 )
beam.SetLoads(       84 ,  0.00018, -0.00011 , -0.00338 )
beam.SetLoads(       85 ,  0.00016, -0.00011 , -0.00325 )
beam.SetLoads(       86 ,  0.00015, -0.00011 , -0.00314 )
beam.SetLoads(       87 ,  0.00014, -0.00010 , -0.00304 )
beam.SetLoads(       88 ,  0.00012, -0.00010 , -0.00295 )
beam.SetLoads(       89 ,  0.00011, -0.00009 , -0.00287 )
beam.SetLoads(       90 ,  0.00010, -0.00007 , -0.00280 )
beam.SetLoads(       91 ,  0.00009, -0.00006 , -0.00274 )
beam.SetLoads(       92 ,  0.00008, -0.00004 , -0.00270 )
beam.SetLoads(       93 ,  0.00008, -0.00002 , -0.00267 )
beam.SetLoads(       94 ,  0.00007,  0.00000 , -0.00265 )
beam.SetLoads(       95 ,  0.00007,  0.00003 , -0.00264 )
beam.SetLoads(       96 ,  0.00006,  0.00006 , -0.00264 )
beam.SetLoads(       97 ,  0.00006,  0.00009 , -0.00265 )
beam.SetLoads(       98 ,  0.00006,  0.00012 , -0.00268 )
beam.SetLoads(       99 ,  0.00006,  0.00016 , -0.00271 )

beam.Run()

beam.ComputeObjectiveFunction(iNode)

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
rigid = 80

for i in range(19,19+rigid):
    node_i = int(beam.elem_py[i].GetNodes()[0, 0] - 1)
    node_j = int(beam.elem_py[i].GetNodes()[1, 0] - 1)
    #plt.plot([ beam.coordinate_X2[node_i],beam.coordinate_X2[node_j] ], [ beam.coordinate_Y2[node_i],beam.coordinate_Y2[node_j] ], [ beam.coordinate_Z2[node_i],beam.coordinate_Z2[node_j] ],'r')
    plt.plot([ beam.coordinate_X1[node_i],beam.coordinate_X1[node_j] ], [ beam.coordinate_Y1[node_i],beam.coordinate_Y1[node_j] ], [ beam.coordinate_Z1[node_i],beam.coordinate_Z1[node_j] ],'b')
    plt.plot([ beam.coordinate_X0[node_i],beam.coordinate_X0[node_j] ], [ beam.coordinate_Y0[node_i],beam.coordinate_Y0[node_j] ], [ beam.coordinate_Z0[node_i],beam.coordinate_Z0[node_j] ],'g')



plt.xlabel('X')
plt.ylabel('Y')
plt.show()

for i in range(beam.nPoint):
    print('{} {}  {}  {} '.format(i, beam.displacement_X[i],beam.displacement_Y[i],beam.displacement_Z[i]) )
