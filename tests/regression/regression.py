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
import numpy as np

beam = CBeamSolver()

loads = [1.0, 2.0, 3.0]
iNode = 20

beam.SetThickness(0.02)
beam.Initialize()
beam.SetLoads(iNode,1,5000)
beam.SetLoads(iNode,2,1000)
beam.Solve()


coordinate_X = []
coordinate_Y = []
coordinate_Z = []

coordinate_X0 = []
coordinate_Y0 = []
coordinate_Z0 = []

for iNode in range(0,21):
    
  coordinate_X.append(beam.ExtractCoordinates(iNode, 0))
  coordinate_Y.append(beam.ExtractCoordinates(iNode, 1))
  coordinate_Z.append(beam.ExtractCoordinates(iNode, 2))  
  
  coordinate_X0.append(beam.ExtractInitialCoordinates(iNode, 0))
  coordinate_Y0.append(beam.ExtractInitialCoordinates(iNode, 1))
  coordinate_Z0.append(beam.ExtractInitialCoordinates(iNode, 2))    

test_val = np.sqrt((coordinate_X[20]-24.020327385028295)**2+
                   (coordinate_Y[20]-16.29552732537537)**2+
                   (coordinate_Z[20]-0.3752371597829022)**2)

print("Tolerance: ",test_val)

# Tolerance is set to 1E-8
if (test_val < 1e-8):
  exit(0)
else:
  exit(1)
