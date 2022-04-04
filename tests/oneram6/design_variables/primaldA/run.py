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

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))

iNode = 19
Loady = -100
Loadz = 7000

########################################
# Initialize and set loads/crossed terms
########################################

primaldA = pyBeamSolver(file_dir, 'configdA.pyBeam')
primaldA.SetLoads(iNode, 0, Loady , Loadz )
primaldA.Run()


posXdA, posYdA, posZdA = primaldA.PrintSolution(iNode)
OFdA = primaldA.ComputeObjectiveFunction(iNode)

exit()
