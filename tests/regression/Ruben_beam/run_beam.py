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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys, os
from pyBeamIO import pyBeamConfig as pyConfig
from pyBeamIO import pyBeamInput as pyInput
import pyBeam

# Load running directory
rundir = os.path.dirname(os.path.realpath(__file__))
confFile = rundir + '/configBeam.cfg'

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
RBE2_py, nRBE2 = pyInput.readRBE2(config['B_MESH']) 
# Parsing Property file
Prop, nProp = pyInput.readProp(config['B_PROPERTY'])

# Initializing objects
beam = pyBeam.CBeamSolver()
inputs = pyBeam.CInput(nPoint, nElem, nRBE2)

# Sending to CInput object 
pyConfig.parseInput(config, inputs, Constr, nConstr)
# Assigning input values to the input object in C++
inputs.SetParameters()
# Initialize the input in the beam solver
beam.InitializeInput(inputs)

# Assigning values to the CNode objects in C++
node = []  
for i in range(nPoint):
   node.append( pyBeam.CNode(node_py[i].GetID()) )
   for j in range(nDim):
      node[i].SetCoordinate(j , float(node_py[i].GetCoord()[j][0]) )
      node[i].SetCoordinate0(j , float(node_py[i].GetCoord0()[j][0]) )
   beam.InitializeNode(node[i], i)
  
# Assigning property values to the property objects in C++
beam_prop = []
for i in range(nProp):
    beam_prop.append(pyBeam.CProperty(i))
    beam_prop[i].SetSectionProperties( Prop[i].GetA(),  Prop[i].GetIyy(),  Prop[i].GetIzz(),  Prop[i].GetJt())
  
# Assigning element values to the element objects in C++ 
element =[]
for i in range(nElem): 
   element.append(pyBeam.CElement(i))
   #element[i].Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input, addouble AuxVector_x, addouble AuxVector_y, addouble AuxVector_z)
   #NB node starts from index 0 and the same happen in beam_prop. But in element_py (connectivity) indexes start from 1 as it is the physical connectivity read from input file
   element[i].Initializer(node[elem_py[i].GetNodes()[0,0] -1], node[elem_py[i].GetNodes()[1,0] -1], beam_prop[elem_py[i].GetProperty() -1], inputs, elem_py[i].GetAuxVector()[0,0], elem_py[i].GetAuxVector()[1,0], elem_py[i].GetAuxVector()[2,0]  )
   beam.InitializeElement(element[i], i)

# IF ANY, assigning RBE2_element values to the RBE2 objects in C++ 
if nRBE2 != 0:
   RBE2 = []
   for i in range(nRBE2): 
      RBE2.append(pyBeam.CRBE2(i))
      RBE2[i].Initializer(node[RBE2_py[i].GetNodes()[0,0] -1], node[RBE2_py[i].GetNodes()[1,0] -1])
      beam.InitializeRBE2(RBE2[i], i)
 
  
  
beam.InitializeStructure()

beam.SetLoads(2 -1,1 - 1,500)
#beam.SetLoads(6,0,100)
#beam.SetLoads(7,0,100)
#beam.SetLoads(8,0,100)

beam.Solve(0)

coordinate_X = []
coordinate_Y = []
coordinate_Z = []

disp_X = []
disp_Y = []
disp_Z = []

coordinate_X0 = []
coordinate_Y0 = []
coordinate_Z0 = []           

for jNode in range(0,2):
    
  coordinate_X.append(beam.ExtractCoordinate(jNode, 0))
  coordinate_Y.append(beam.ExtractCoordinate(jNode, 1))
  coordinate_Z.append(beam.ExtractCoordinate(jNode, 2))
  
  disp_X.append(beam.ExtractDisplacements(jNode, 0))
  disp_Y.append(beam.ExtractDisplacements(jNode, 1))
  disp_Z.append(beam.ExtractDisplacements(jNode, 2))  
  
  coordinate_X0.append(beam.ExtractCoordinate0(jNode, 0))
  coordinate_Y0.append(beam.ExtractCoordinate0(jNode, 1))
  coordinate_Z0.append(beam.ExtractCoordinate0(jNode, 2))
  
test_val = np.sqrt((disp_X[1])**2+
                   (disp_Y[1])**2+
                   (disp_Z[1])**2)
                   
print("Displacement tip: ",test_val)          
  
fig = plt.figure()
ax = fig.add_subplot(111)
  
plt.plot(coordinate_X, coordinate_Y)
plt.plot(coordinate_X0, coordinate_Y0)
plt.xlabel('X')
plt.ylabel('Y')
plt.xlim((-0.01, 0.01)) 
plt.show()

