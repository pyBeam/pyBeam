#!/usr/bin/env python
#
# pyBeam, a Beam Solver
#
# Copyright (C) 2018 Rocco Bombardieri, Ruben Sanchez, Rauno Cavallaro
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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys, os
from pyBeamLib import pyBeamSolver

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))

beam = pyBeamSolver(file_dir ,'config_NL.cfg')

beam.SetLoads(20,0,0,50000)

beam.Run()

beam.PrintDisplacements(20)

beam.Restart()

beam.PrintDisplacements(20)

success = beam.TestNodePosition(20,1e-8,19.05326647582835,12.150042136520977,16.32312936016357)



coordinate_X = []
coordinate_Y = []
coordinate_Z = []
coordinate_X0 = []
coordinate_Y0 = []
coordinate_Z0 = []


for jNode in range(0,beam.nPoint):

  dispX, dispY, dispZ = beam.ExtractDisplacements(jNode)

  coordX, coordY, coordZ = beam.GetInitialCoordinates(jNode)

  coordinate_X.append(coordX + dispX)
  coordinate_Y.append(coordY + dispY)
  coordinate_Z.append(coordZ + dispZ)
  
  coordinate_X0.append(coordX)
  coordinate_Y0.append(coordY)
  coordinate_Z0.append(coordZ)
  

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
   

  
plt.plot(coordinate_X, coordinate_Y, coordinate_Z)
plt.plot(coordinate_X0, coordinate_Y0, coordinate_Z0)
plt.xlabel('X')
plt.ylabel('Y')
plt.show()




exit(success)
