#!/usr/bin/env python3
#
# pyBeam, an open-source Beam Solver
#
# Copyright (C) 2019 by the authors
# 
# File developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
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

import shutil
import numpy as np
import sys, os

from pyBeamLib import pyBeamSolver
from pyBeamLibAD import pyBeamSolverAD

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))

# Set test node
iNode = 20

############################
# Run Primal Case
############################
primal = pyBeamSolver(file_dir, 'config.pyBeam')
primal.SetLoads(iNode, 0, -1000, 20000)
primal.Run()
posX, posY, posZ = primal.PrintSolution(iNode)
primal.ComputeObjectiveFunction(iNode)

# Move restart file into solution file before running FD
shutil.move("restart.pyBeam", "solution.pyBeam")

############################
# Compute Finite Differences
############################
primal.SetLowVerbosity()

# Load X
increment = 0.1
primal.SetLoads(iNode, 0+increment, -1000, 20000)
primal.Run()
of_xpos = primal.ComputeObjectiveFunction(iNode)
primal.SetLoads(iNode, 0-increment, -1000, 20000)
primal.Run()
of_xneg = primal.ComputeObjectiveFunction(iNode)
sensX_FD = (of_xpos - of_xneg) / (2*increment)

# Load Y
increment = 1
primal.SetLoads(iNode, 0, -1000+increment, 20000)
primal.Run()
of_ypos = primal.ComputeObjectiveFunction(iNode)
primal.SetLoads(iNode, 0, -1000-increment, 20000)
primal.Run()
of_yneg = primal.ComputeObjectiveFunction(iNode)
sensY_FD = (of_ypos - of_yneg) / (2*increment)

# Load Z
increment = 10
primal.SetLoads(iNode, 0, -1000, 20000+increment)
primal.Run()
of_zpos = primal.ComputeObjectiveFunction(iNode)
primal.SetLoads(iNode, 0, -1000, 20000-increment)
primal.Run()
of_zneg = primal.ComputeObjectiveFunction(iNode)
sensZ_FD = (of_zpos - of_zneg) / (2*increment)

############################
# Run Adjoint Case
############################
adjoint = pyBeamSolverAD(file_dir, 'configAD.pyBeam')
adjoint.SetLoads(iNode, 0, -1000, 20000)
adjoint.ReadRestart() 
adjoint.StartRecording()
adjoint.SetDependencies()
adjoint.Restart()
adjoint.ComputeObjectiveFunction(iNode)
adjoint.StopRecording()
adjoint.ComputeAdjoint()

sensE = adjoint.PrintSensitivityE()
sensX, sensY, sensZ = adjoint.PrintSensitivityLoad(iNode)

############################
# Tests
############################
print("\n\n############################\n TESTS \n############################\n")
test_nodePos = primal.TestNodePosition(posX, posY, posZ, 28.177002149725, -6.121758821772, 7.000638871880, 1E-8)
test_sensLoad = adjoint.TestLoadSensitivity(iNode, sensX, sensY, sensZ, sensX_FD, sensY_FD, sensZ_FD)
test_sensE = adjoint.TestSensitivityE(sensE, -1.6748137651862045e-10)

test = test_nodePos + test_sensLoad + test_sensE
exit(test)


