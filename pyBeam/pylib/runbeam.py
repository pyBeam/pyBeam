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
import sys
sys.path.append('..')
import in_out
import swig
#from pyBeamAD import CBeamSolver 

def RunBeam():

  confFile = '/media/rocco/290CF0EB732D1122/Adjoint_project/pyBeam/tests/OneraM6/BEAM_config.cfg'

# Parsing Conf file
BEAM_config = in_out.BEAMConfig(confFile) 		# FSI configuration file

# Parsing mesh file
nDim = in_out.readDimension(BEAM_config['B_MESH'])
node_py, nPoint = in_out.readMesh(BEAM_config['B_MESH'],nDim)
elem_py, nElem = in_out.readConnectivity( BEAM_config['B_MESH']) 
Constr, nConstr = in_out.readConstr(BEAM_config['B_MESH'])  
# Parsing Property file
Prop, nProp = in_out.readProp(BEAM_config['B_PROPERTY'])


# Initializing objects
beam = swig.CBeamSolver()
inputs = swig.CInput()
  
    
# Sending to CInput object 
swig.Input_parsing(BEAM_config, inputs, Constr, nConstr)
# Assigning input values to the input object in C++
inputs.SetParameters()

# Assigning values to the CNode objects in C++
node = []  
for i in range(nPoint):
   node.append( swig.CNode(node_py[i].GetID()) )
   for j in range(nDim):
      node[i].SetCoordinate(j , float(node_py[i].GetCoord()[j][0]) )
      node[i].SetCoordinate0(j , node_py[i].GetCoord0()[j][0])
    
# Assigning property values to the property objects in C++
beam_prop = []
for i in range(nProp):
    beam_prop.append(swig.CProperty(i))
    beam_prop[i].SetSectionProperties( Prop[i].GetA(),  Prop[i].GetIyy(),  Prop[i].GetIzz(),  Prop[i].GetJt())
  
# Assigning element values to the property objects in C++ 
element =[]
for i in range(nElem): 
   element.append(swig.CElement(i))
   #element[i].Initializer(CNode* Node1, CNode* Node2, CProperty* Property, CInput* Input, addouble AuxVector_x, addouble AuxVector_y, addouble AuxVector_z)
   #NB node starts from index 0 and the same happen in beam_prop. But in element_py (connectivity) indexes start from 1 as it is the physical connectivity read from input file
   element[i].Initializer(node[elem_py[i].GetNodes()[0,0] -1], node[elem_py[i].GetNodes()[1,0] -1], beam_prop[elem_py[i].GetProperty() -1], inputs, elem_py[i].GetAuxVector()[0,0], elem_py[i].GetAuxVector()[1,0], elem_py[i].GetAuxVector()[2,0]  )
   
  
  iNode = 20

  #beam.SetThickness(0.02)
  beam.StartRecording()
  beam.Initialize(inputs)
  
  beam.SetLoads(iNode,1,5000)
  beam.SetLoads(iNode,2,1000)
  beam.RegisterLoads()
  beam.Solve()
  displacement = beam.OF_NodeDisplacement(iNode)
  beam.StopRecording()

  thickness_gradient = beam.ComputeAdjoint()

  print("Objective Function - Displacement(",iNode,") = ", displacement)
  print("t' = ", thickness_gradient)

  coordinate_X = []
  coordinate_Y = []
  coordinate_Z = []

  coordinate_X0 = []
  coordinate_Y0 = []
  coordinate_Z0 = []

  for iNode in range(0,21):
    
    print("F'(",iNode,") = (", beam.ExtractLoadGradient(iNode,0), beam.ExtractLoadGradient(iNode,1), beam.ExtractLoadGradient(iNode,2), ")")
    
    coordinate_X.append(beam.ExtractCoordinates(iNode, 0))
    coordinate_Y.append(beam.ExtractCoordinates(iNode, 1))
    coordinate_Z.append(beam.ExtractCoordinates(iNode, 2))  
    
    coordinate_X0.append(beam.ExtractInitialCoordinates(iNode, 0))
    coordinate_Y0.append(beam.ExtractInitialCoordinates(iNode, 1))
    coordinate_Z0.append(beam.ExtractInitialCoordinates(iNode, 2))    
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  
  plt.plot(coordinate_X, coordinate_Y, coordinate_Z)
  plt.plot(coordinate_X0, coordinate_Y0, coordinate_Z0)
  plt.show()

  
  
def Input_parsing(BEAM_config, inputs,Constr, nConstr):  
    
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
    
    #Now setting the constraints
    inputs.SetnConstr(nConstr)
    for iConstr in range(nConstr): 
        inputs.SetSingleConstr( iConstr, Constr[iConstr,0], Constr[iConstr,1] )
  
# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    RunBeam()                   
                   
                   

                   
                   
                   
                   
                     
