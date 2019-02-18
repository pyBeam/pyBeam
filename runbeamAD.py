#
# pyBeam, a Beam Solver
#
# Copyright (C) 2018 Ruben Sanchez, Rauno Cavallaro
# 
# Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
#             Rauno Cavallaro (Carlos III University Madrid)
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


from pyBeamAD import CBeamSolver 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

beam = CBeamSolver()

loads = [1.0, 2.0, 3.0]
iNode = 20

beam.SetThickness(0.02)
beam.StartRecording()
beam.Initialize()
beam.SetLoads(iNode,1,5000)
beam.SetLoads(iNode,2,1000)
beam.RegisterLoads()
beam.Solve()
displacement = beam.OF_NodeDisplacement(iNode)
beam.StopRecording()

thickness_gradient = beam.ComputeAdjoint()

print("Objective Function - Displacement(",iNode,") = ", displacement)
print("t' = ", thickness_gradient)

coordinate_X = []
coordinate_Y = []
coordinate_Z = []

coordinate_X0 = []
coordinate_Y0 = []
coordinate_Z0 = []

for iNode in range(0,21):
  
  print("F'(",iNode,") = (", beam.ExtractLoadGradient(iNode,0), beam.ExtractLoadGradient(iNode,1), beam.ExtractLoadGradient(iNode,2), ")")
  
  coordinate_X.append(beam.ExtractCoordinates(iNode, 0))
  coordinate_Y.append(beam.ExtractCoordinates(iNode, 1))
  coordinate_Z.append(beam.ExtractCoordinates(iNode, 2))  
  
  coordinate_X0.append(beam.ExtractInitialCoordinates(iNode, 0))
  coordinate_Y0.append(beam.ExtractInitialCoordinates(iNode, 1))
  coordinate_Z0.append(beam.ExtractInitialCoordinates(iNode, 2))    

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plt.plot(coordinate_X, coordinate_Y, coordinate_Z)
plt.plot(coordinate_X0, coordinate_Y0, coordinate_Z0)
plt.show(block=False)

wait = input("Press Enter to finalize.")
