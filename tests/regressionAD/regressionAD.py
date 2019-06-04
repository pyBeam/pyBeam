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
from pyBeamLib import pyBeamSolver
import pyBeamAD

# Load running directory
rundir = os.path.dirname(os.path.realpath(__file__))
confFile = rundir + '/config.cfg'

# Parsing Conf file
config = pyConfig.pyBeamConfig(confFile)  # Beam configuration file

# Specifically added for the test
config['B_PROPERTY'] = rundir + '/' + config['B_PROPERTY'][:]
config['B_MESH'] = rundir + '/' + config['B_MESH'][:]
  
# Parsing mesh file
nDim = pyInput.readDimension(config['B_MESH'])
node_py, nPoint = pyInput.readMesh(config['B_MESH'],nDim)
elem_py, nElem = pyInput.readConnectivity( config['B_MESH'])
Constr, nConstr = pyInput.readConstr(config['B_MESH'])
# Parsing Property file
Prop, nProp = pyInput.readProp(config['B_PROPERTY'])

# Initializing objects
beam = pyBeamAD.CBeamSolver()
inputs = pyBeamAD.CInput(nPoint, nElem)

# Start recording
beam.StartRecording()
    
# Sending to CInput object 
pyConfig.parseInput(config, inputs, Constr, nConstr)
# Assigning input values to the input object in C++
inputs.SetParameters()
# Initialize the input in the beam solver
beam.InitializeInput(inputs)

# Assigning values to the CNode objects in C++
node = []  
for i in range(nPoint):
   node.append( pyBeamAD.CNode(node_py[i].GetID()) )
   for j in range(nDim):
      node[i].SetCoordinate(j , float(node_py[i].GetCoord()[j][0]) )
      node[i].SetCoordinate0(j , float(node_py[i].GetCoord0()[j][0]) )
   beam.InitializeNode(node[i], i)
    
# Assigning property values to the property objects in C++
beam_prop = []
for i in range(nProp):
    beam_prop.append(pyBeamAD.CProperty(i))
    beam_prop[i].SetSectionProperties( Prop[i].GetA(),  Prop[i].GetIyy(),  Prop[i].GetIzz(),  Prop[i].GetJt())
  
# Assigning element values to the property objects in C++ 
element =[]
for i in range(nElem): 
   element.append(pyBeamAD.CElement(i))
   #element[i].Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input, addouble AuxVector_x, addouble AuxVector_y, addouble AuxVector_z)
   #NB node starts from index 0 and the same happen in beam_prop. But in element_py (connectivity) indexes start from 1 as it is the physical connectivity read from input file
   element[i].Initializer(node[elem_py[i].GetNodes()[0,0] -1], node[elem_py[i].GetNodes()[1,0] -1], beam_prop[elem_py[i].GetProperty() -1], inputs, elem_py[i].GetAuxVector()[0,0], elem_py[i].GetAuxVector()[1,0], elem_py[i].GetAuxVector()[2,0]  )
   beam.InitializeElement(element[i], i)
  
beam.InitializeStructure()
  
iNode = 20
beam.SetLoads(iNode,1,5000)
beam.SetLoads(iNode,2,1000)
beam.RegisterLoads()
beam.Solve(0)
displacement = beam.OF_NodeDisplacement(iNode)
beam.StopRecording()

beam.ComputeAdjoint()

print("Objective Function - Displacement(",iNode,") = ", displacement)

coordinate_X = []
coordinate_Y = []
coordinate_Z = []

coordinate_X0 = []
coordinate_Y0 = []
coordinate_Z0 = []

for iNode in range(0,21):
  
  print("F'(",iNode,") = (", beam.ExtractLoadGradient(iNode,0), beam.ExtractLoadGradient(iNode,1), beam.ExtractLoadGradient(iNode,2), ")")
  
#  coordinate_X.append(beam.ExtractCoordinates(iNode, 0))
#  coordinate_Y.append(beam.ExtractCoordinates(iNode, 1))
#  coordinate_Z.append(beam.ExtractCoordinates(iNode, 2))  
  
#  coordinate_X0.append(beam.ExtractInitialCoordinates(iNode, 0))
#  coordinate_Y0.append(beam.ExtractInitialCoordinates(iNode, 1))
#  coordinate_Z0.append(beam.ExtractInitialCoordinates(iNode, 2))    

delta_gradient_x = beam.ExtractLoadGradient(iNode,0) + 0.0017035445423928379
delta_gradient_y = beam.ExtractLoadGradient(iNode,1) - 0.0020190915772016127
delta_gradient_z = beam.ExtractLoadGradient(iNode,2) - 1.4059868175534144e-05

test_val = np.sqrt(delta_gradient_x**2 + delta_gradient_y**2 + delta_gradient_z**2)

print("Tolerance: ",test_val)


# Tolerance is set to 1E-8
if (abs(test_val) < abs(1e-8)):
  exit(0)
else:
  exit(1)
