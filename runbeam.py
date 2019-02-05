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


from pyBeam import CBeamSolver 
import matplotlib.pyplot as plt

testobject = CBeamSolver()

loads = [1.0, 2.0, 3.0]

testobject.Initialize()
testobject.SetLoads(20,1,5000)
testobject.Solve()

coordinate_X = []
coordinate_Y = []

coordinate_X0 = []
coordinate_Y0 = []

for iNode in range(0,21):
  
  coordinate_X.append(testobject.ExtractCoordinates(iNode, 0))
  coordinate_Y.append(testobject.ExtractCoordinates(iNode, 1))
  
  coordinate_X0.append(testobject.ExtractInitialCoordinates(iNode, 0))
  coordinate_Y0.append(testobject.ExtractInitialCoordinates(iNode, 1))

plt.plot(coordinate_X, coordinate_Y)
plt.plot(coordinate_X0, coordinate_Y0)
plt.show()
