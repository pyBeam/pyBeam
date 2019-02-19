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


import pdb
import numpy as np
import sys, os
sys.path.append(str(os.path.realpath(__file__))[:-31] + '/pyBeam')
sys.path.append(str(os.path.realpath(__file__))[:-31] + '/pyBeam/python')
#sys.path.append('../../pyBeam')
#sys.path.append('../../pyBeam/python')
import in_out
import swig

def Input_parsing(BEAM_config, inputs):  
    
    inputs.SetBeamLength(BEAM_config['B_LENGTH'])
    inputs.SetWebThickness(BEAM_config['W_THICKNESS'])
    inputs.SetWebHeight(BEAM_config['W_HEIGHT'])
    inputs.SetFlangeWidth(BEAM_config['F_WIDTH'])
    inputs.SetYoungModulus(BEAM_config['Y_MODULUS'])
    inputs.SetPoisson(BEAM_config['POISSON'])
    inputs.SetDensity(BEAM_config['RHO'])
    inputs.SetLoad(BEAM_config['LOAD'])
    inputs.SetFollowerFlag(BEAM_config['FOLLOWER_FLAG'])
    inputs.SetLoadSteps(BEAM_config['LOAD_STEPS'])
    inputs.SetNStructIter(BEAM_config['N_STRUCT_ITER'])
    inputs.SetConvCriterium(BEAM_config['CONV_CRITERIUM'])  



confFile = str(os.path.realpath(__file__))[:-25] + '/OneraM6/BEAM_config.cfg'
BEAM_config = in_out.BEAMConfig(confFile) 		# FSI configuration file

# Initializing objects
beam = swig.CBeamSolver()
inputs = swig.CInput()
  
    
# Parsing config file ans sending to CInput object  
Input_parsing(BEAM_config, inputs)
inputs.SetParameters()
  
# Specifically added for the test
BEAM_config['B_PROPERTY'] = str(os.path.realpath(__file__))[:-25] + BEAM_config['B_PROPERTY'][2:]
BEAM_config['B_MESH'] = str(os.path.realpath(__file__))[:-25] + BEAM_config['B_MESH'][2:]  
# Parsing mesh file
nDim = in_out.readDimension(BEAM_config['B_MESH'])
node, nPoint = in_out.readMesh(BEAM_config['B_MESH'],nDim)
Elem, nElem = in_out.readConnectivity( BEAM_config['B_MESH'])  
# Parsing Property file
Prop, nProp = in_out.readProp(BEAM_config['B_PROPERTY'])

# Assigning property values to the property objects in C++
beam_prop = []
for i in range(nProp):
    beam_prop.append(swig.CProperty(i))
    beam_prop[i].SetSectionProperties( Prop[i].GetA(),  Prop[i].GetIyy(),  Prop[i].GetIzz(),  Prop[i].GetJt())
 
iNode = 20

beam.Initialize(inputs)
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
  
  
  
  
  
  
