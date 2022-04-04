#!/usr/bin/env python
#
# pyBeam, an open-source Beam Solver
#
# Copyright (C) 2019 by the authors
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


sys.path.append('/home/rauno/pyBeamSU2/bin')

from pyBeamLib import pyBeamSolver
from pyBeamLibAD import pyBeamSolverAD

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))

iNode = 19
Loady = -100
Loadz = 7000

########################################
# Initialize and set loads/crossed terms
########################################
#os.chdir('./nominal')
primal = pyBeamSolver(file_dir, 'config.pyBeam')
primal.SetLoads(iNode, 0, Loady , Loadz )
primal.Run()

# Copy the solution file to the running directory
#solution_file = file_dir + '/restart.pyBeam'
#shutil.move(solution_file, "solution.pyBeam")

posX, posY, posZ = primal.PrintSolution(iNode)
OF = primal.ComputeObjectiveFunction(iNode)
exit()

###  dA
primaldA = pyBeamSolver(file_dir, 'configdA.pyBeam')

iNode = 19
primaldA.SetLoads(iNode, 0, -100, 7000)
  
primaldA.Run()


posX, posY, posZ = primaldA.PrintSolution(iNode)
OFdA = primaldA.ComputeObjectiveFunction(iNode)

deriv = (OFdA-OF)/(1.0E-7)
print("sensitivity FD = ",deriv)


########################################
# Initialize and set loads/crossed terms
########################################

adjoint = pyBeamSolverAD(file_dir, 'configAD.pyBeam')


adjoint.SetLoads(iNode, 0, -100, 7000)


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
print("4-Eval objective!!")
adjoint.ComputeObjectiveFunction(iNode)
print("5-Stopping recording!!")
adjoint.StopRecording()
print("6-COmputing Adjoint")
adjoint.ComputeAdjoint()

adjoint.PrintSensitivityDV()
sensE = adjoint.PrintSensitivityE()
adjoint.PrintSensitivityLoad(iNode)

exit()
