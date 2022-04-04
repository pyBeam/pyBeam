#!/usr/bin/env python
#
# pyBeam, an open-source Beam Solver
#
# Copyright (C) 2019 by the authors
# 
# File Developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
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
import sys, os, csv

#from pyBeamLib import pyBeamSolver
from pyBeamLibAD import pyBeamSolverAD

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))


########################################
# Initialize and set loads/crossed terms
########################################

adjoint = pyBeamSolverAD(file_dir, 'configAD.pyBeam')

iNode = 99
adjoint.SetLoads(iNode, 0, -1000, 20000)
#  adjoint.SetDisplacementAdjoint(iNode, adjX[iNode], adjY[iNode], adjZ[iNode])
 

############################
# Solve adjoint
############################
  
adjoint.ReadRestart() 
print("1-Start Recording!!")
adjoint.StartRecording()
print("Done!")
print("2-Setting Dependencies!!")
adjoint.SetDependencies()
print("3-Restarting 1 ieration to record dependencies!!")
adjoint.Restart()
print("4-Stopping recording!!")
adjoint.StopRecording()
print("5-COmputing Adjoint")
adjoint.ComputeAdjoint()

#===
sensE = adjoint.PrintSensitivityE()
adjoint.PrintSensitivityLoad(99)
adjoint.PrintSensitivityDV()
 
exit()
############################
# Tests
############################

print("\n############################\n TEST \n############################\n")
test = adjoint.TestSensitivityE(sensE, 1.0546216892681002e-13)

exit(test)

