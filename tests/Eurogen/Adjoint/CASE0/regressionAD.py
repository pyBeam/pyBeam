#!/usr/bin/env python
#
# pyBeam, a Beam Solver
#
# Copyright (C) 2019 Rocco Bombardieri, Ruben Sanchez , Rauno Cavallaro
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


import sys, os
from pyBeamLibAD import pyBeamSolverAD

# Load running directory
file_dir = os.path.dirname(os.path.realpath(__file__))

beam = pyBeamSolverAD(file_dir, 'config.cfg')

iNode = 21 - 1

beam.SetLoads(iNode,0,0,60000)

beam.StartRecording()

beam.SetDependencies()

#beam.SetDisplacementAdjoint(iNode, 0.0, 0.0, 0.0)

beam.Run() 

beam.ComputeObjectiveFunction(iNode)

beam.StopRecording()

beam.ComputeAdjoint()

beam.PrintSensitivitiesAllLoads()

