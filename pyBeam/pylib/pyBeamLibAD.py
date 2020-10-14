#!/usr/bin/env python
#
# pyBeam, a Beam Solver
#
# Copyright (C) 2019 by the authors
# 
# File developers: Rocco Bombardieri (Carlos III University Madrid)
#                  Ruben Sanchez (SciComp, TU Kaiserslautern)
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# Imports from pyBeam
from pyBeamIO import pyBeamConfig as pyConfig
from pyBeamIO import pyBeamInput as pyInput
import numpy as np
import pyBeamAD

# ----------------------------------------------------------------------
#  Beam object
# ----------------------------------------------------------------------

class pyBeamSolverAD:
  """Description"""

  def __init__(self, file_dir, config_fileName):
    """ Description. """

    self.file_dir = file_dir
    self.Config_file = self.file_dir + '/' + config_fileName
    self.Config = {}

    print("\n---------------------------------------------------------------------------")
    print("|                                                                         |")
    print("| pyBeam, a Beam Solver (AD) - Release 0.1 (beta)                         |")
    print("|                            - https://github.com/rsanfer/pyBeam          |")
    print("|                                                                         |")
    print("---------------------------------------------------------------------------")
    print("|                                                                         |")
    print("| Copyright 2018-2019 Rocco Bombardieri (Carlos III University Madrid)    |")
    print("|                     Rauno Cavallaro (Carlos III University Madrid)      |")
    print("|                     Ruben Sanchez (SciComp, TU Kaiserslautern)          |")
    print("|                     Tim Albring (SciComp, TU Kaiserslautern)            |")
    print("|                                                                         |")
    print("| pyBeam is free software: you can redistribute it and/or                 |")
    print("| modify it under the terms of the GNU Affero General Public License      |")
    print("| as published by the Free Software Foundation, either version 3 of the   |")
    print("| License, or (at your option) any later version.                         |")
    print("|                                                                         |")
    print("| pyBeam is distributed in the hope that it will be useful,               |")
    print("| but WITHOUT ANY WARRANTY; without even the implied warranty             |")
    print("| of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                 |")
    print("| See the GNU Affero General Public License for more details.             |")
    print("|                                                                         |")
    print("| You should have received a copy of the GNU Affero                       |")
    print("| General Public License along with pyBeam.                               |")
    print("| If not, see <http://www.gnu.org/licenses/>.                             |")
    print("|                                                                         |")
    print("---------------------------------------------------------------------------\n")

    # Parsing config file
    self.Config = pyConfig.pyBeamConfig(self.Config_file)  # Beam configuration file

    self.Mesh_file = self.file_dir + '/' + self.Config['MESH_FILE']
    self.Property = self.file_dir + '/' + self.Config['PROPERTY_FILE']

    # Parsing mesh file
    self.nDim= pyInput.readDimension(self.Mesh_file)
    self.node_py, self.nPoint = pyInput.readMesh(self.Mesh_file, self.nDim)
    self.elem_py, self.nElem = pyInput.readConnectivity(self.Mesh_file)
    self.Constr, self.nConstr = pyInput.readConstr(self.Mesh_file)
    self.RBE2_py, self.nRBE2 = pyInput.readRBE2(self.Mesh_file)
    # Parsing Property file
    self.Prop, self.nProp = pyInput.readProp(self.Property)

    # Initializing objects
    self.beam = pyBeamAD.CBeamSolver()
    self.inputs = pyBeamAD.CInput(self.nPoint, self.nElem, self.nRBE2, self.nProp)

    # Start recording
    print("--> Initialization successful!")

    # Sending to CInput object
    pyConfig.parseInput(self.Config, self.inputs, self.Constr, self.nConstr)
    # Assigning input values to the input object in C++
    self.inputs.SetParameters()
    # Set the discrete adjoint flag to true
    self.inputs.SetDiscreteAdjoint()
    # Initialize the input in the beam solver
    self.beam.InitializeInput(self.inputs)

    # Assigning values to the CNode objects in C++
    self.node = []
    for i in range(self.nPoint):
        self.node.append(pyBeamAD.CNode(self.node_py[i].GetID()))
        for j in range(self.nDim):
            self.node[i].InitCoordinate(j, float(self.node_py[i].GetCoord()[j][0]))
        self.beam.InitializeNode(self.node[i], i)

    # Assigning property values to the property objects in C++
    self.beam_prop = []
    self.nPropDVs = 0
    for i in range(self.nProp):
        self.beam_prop.append(pyBeamAD.CProperty(i))
        if self.Prop[i].GetFormat() == "N":
            self.nPropDVs = self.nPropDVs + 4
            self.beam_prop[i].SetSectionProperties(self.Prop[i].GetA(), self.Prop[i].GetIyy(), self.Prop[i].GetIzz(), self.Prop[i].GetJt())
        elif self.Prop[i].GetFormat() == "S":
            self.nPropDVs = self.nPropDVs + 6
            self.beam_prop[i].SetSectionProperties(self.Prop[i].GetC_wb(), self.Prop[i].Geth(), self.Prop[i].Gett_sk(), self.Prop[i].Gett_sp(),\
            self.Prop[i].GetA_fl(), self.Prop[i].Getn_stiff(),self.Prop[i].GetA_stiff() )
            #print(self.Prop[i].GetA())
        else:
            raise ValueError("Unknown paramter for Property CARD input. Execution aborted")
        self.beam.InitializeProp(self.beam_prop[i], i)

    # Assigning element values to the element objects in C++
    self.element = []
    for i in range(self.nElem):
        self.element.append(pyBeamAD.CElement(i))
        self.element[i].Initializer(self.node[self.elem_py[i].GetNodes()[0, 0] - 1], self.node[self.elem_py[i].GetNodes()[1, 0] - 1],
                               self.beam_prop[self.elem_py[i].GetProperty() - 1], self.inputs, self.elem_py[i].GetAuxVector()[0, 0],
                               self.elem_py[i].GetAuxVector()[1, 0], self.elem_py[i].GetAuxVector()[2, 0])
        self.beam.InitializeElement(self.element[i], i)

    # Here we need to pass the AeroPoint matrix of the wing grid
    # IF ANY, assigning RBE2_element values to the RBE2 objects in C++
    if self.nRBE2 != 0:
        self.RBE2 = []
        for i in range(self.nRBE2):
            self.RBE2.append(pyBeamAD.CRBE2(i))
            self.RBE2[i].Initializer(self.node[self.RBE2_py[i].GetNodes()[0, 0] - 1], self.node[self.RBE2_py[i].GetNodes()[1, 0] - 1])
            self.beam.InitializeRBE2(self.RBE2[i], i)


    # Initialize structures to store the coordinates and displacements

    self.coordinate_X = []
    self.coordinate_Y = []
    self.coordinate_Z = []

    self.coordinate_X0 = []
    self.coordinate_Y0 = []
    self.coordinate_Z0 = []

    self.displacement_X = []
    self.displacement_Y = []
    self.displacement_Z = []

    # finally intializing the structure for the solver
    self.beam.InitializeStructure()

    print("--> Initialization successful")
    print("\n---------------------------------------------------------------------------\n")


  def RegisterLoads(self):
    """ This function starts load registration for AD  """
    self.beam.RegisterLoads()

  def StartRecording(self):
    """ This function stops registration for AD  """
    self.beam.StartRecording()

  def SetDependencies(self):
    """ This function stops registration for AD  """
    self.beam.SetDependencies()

  def StopRecording(self):
    """ This function stops registration for AD  """
    self.beam.StopRecording()
    
  def StopRecordingWeight(self):
    """ This function stops registration for AD  """
    self.beam.StopRecordingWeight()
    
  def StopRecordingKS(self):
    """ This function stops registration for AD  """
    self.beam.StopRecordingKS()
    
  def StopRecordingEA(self):
    """ This function stops registration for AD  """
    self.beam.StopRecordingEA()
    
  def StopRecordingNint(self):
    """ This function stops registration for AD  """
    self.beam.StopRecordingNint()
    
    
  def ComputeAdjoint(self):
    """ This function computes Adjoint for AD  """
    self.beam.ComputeAdjoint()
    
  def ComputeAdjointNint(self):
    """ This function computes Adjoint for AD  """
    self.beam.ComputeAdjointNint()
    
  def ComputeAdjointWeight(self):
    """ This function computes Adjoint for AD  """
    self.beam.ComputeAdjointWeight()
    
  def ComputeAdjointKS(self):
    """ This function computes Adjoint for AD  """
    self.beam.ComputeAdjointKS()
    
  def ComputeAdjointEA(self):
    """ This function computes Adjoint for AD  """
    self.beam.ComputeAdjointEA()
    
    
    
  def SetLoads(self, iVertex, loadX, loadY, loadZ):

    """ This function sets the load  """
    self.beam.SetLoads(iVertex, 0, loadX)
    self.beam.SetLoads(iVertex, 1, loadY)
    self.beam.SetLoads(iVertex, 2, loadZ)

  def GetLoadSensitivity(self, iVertex):

    """ This function returns the load sensitivity  """
    sensX = self.beam.ExtractLoadGradient(iVertex, 0)
    sensY = self.beam.ExtractLoadGradient(iVertex, 1)
    sensZ = self.beam.ExtractLoadGradient(iVertex, 2)

    return sensX, sensY, sensZ

  def SetDisplacementAdjoint(self, iVertex, adjX, adjY, adjZ):
    """ This function sets the load  """
    self.beam.StoreDisplacementAdjoint(iVertex, 0, adjX)
    self.beam.StoreDisplacementAdjoint(iVertex, 1, adjY)
    self.beam.StoreDisplacementAdjoint(iVertex, 2, adjZ)


  def ComputeObjectiveFunction(self, iNode):

    """ This function computes the objective function (Important to be recorded) """
    displacement = self.beam.OF_NodeDisplacement(iNode)
    print("Objective Function - Displacement(", iNode, ") = ", displacement)
    return displacement

  
  #####  DEBUG
  def ComputeNint(self):

    Nint = self.beam.EvalNint()
    return Nint


  def ComputeWeight(self):
    """ This function computes the response weight of the structure (important to be recorded) """
    weight = self.beam.EvalWeight()
    return weight


  def ComputeResponseKSStress(self):
    """ This function computes the KS stress on the structure (important to be recorded) """
    KSStress= self.beam.EvalKSStress()
    return KSStress

  def ComputeEA(self):
    EA = self.beam.EvalEA()
    return EA


  def Run(self):
    """ This function runs the solver and stores the results.
        Needs to be run after __SetLoads """

    self.beam.Solve(0)

    self.coordinate_X = []
    self.coordinate_Y = []
    self.coordinate_Z = []

    self.coordinate_X0 = []
    self.coordinate_Y0 = []
    self.coordinate_Z0 = []

    self.displacement_X = []
    self.displacement_Y = []
    self.displacement_Z = []

    for jNode in range(0,self.nPoint):

        self.coordinate_X.append(self.beam.ExtractCoordinate(jNode, 0))
        self.coordinate_Y.append(self.beam.ExtractCoordinate(jNode, 1))
        self.coordinate_Z.append(self.beam.ExtractCoordinate(jNode, 2))

        self.coordinate_X0.append(self.beam.ExtractCoordinate0(jNode, 0))
        self.coordinate_Y0.append(self.beam.ExtractCoordinate0(jNode, 1))
        self.coordinate_Z0.append(self.beam.ExtractCoordinate0(jNode, 2))

        self.displacement_X.append(self.beam.ExtractDisplacements(jNode, 0))
        self.displacement_Y.append(self.beam.ExtractDisplacements(jNode, 1))
        self.displacement_Z.append(self.beam.ExtractDisplacements(jNode, 2))


  def RunLin(self):
    """ This function runs the solver and stores the results.
        Needs to be run after __SetLoads """

    self.beam.SolveLin(0)

    self.coordinate_X = []
    self.coordinate_Y = []
    self.coordinate_Z = []

    self.coordinate_X0 = []
    self.coordinate_Y0 = []
    self.coordinate_Z0 = []

    self.displacement_X = []
    self.displacement_Y = []
    self.displacement_Z = []

    for jNode in range(0,self.nPoint):

        self.coordinate_X.append(self.beam.ExtractCoordinate(jNode, 0))
        self.coordinate_Y.append(self.beam.ExtractCoordinate(jNode, 1))
        self.coordinate_Z.append(self.beam.ExtractCoordinate(jNode, 2))

        self.coordinate_X0.append(self.beam.ExtractCoordinate0(jNode, 0))
        self.coordinate_Y0.append(self.beam.ExtractCoordinate0(jNode, 1))
        self.coordinate_Z0.append(self.beam.ExtractCoordinate0(jNode, 2))

        self.displacement_X.append(self.beam.ExtractDisplacements(jNode, 0))
        self.displacement_Y.append(self.beam.ExtractDisplacements(jNode, 1))
        self.displacement_Z.append(self.beam.ExtractDisplacements(jNode, 2))
        
        
  def ReadRestart(self):

      self.beam.ReadRestart()

  def Restart(self):
      """ This function runs the restart and stores the results.
          Needs to be run after __SetLoads """

      self.beam.RunRestart(0)

      self.coordinate_X = []
      self.coordinate_Y = []
      self.coordinate_Z = []

      self.coordinate_X0 = []
      self.coordinate_Y0 = []
      self.coordinate_Z0 = []

      self.displacement_X = []
      self.displacement_Y = []
      self.displacement_Z = []

      for jNode in range(0, self.nPoint):
          self.coordinate_X.append(self.beam.ExtractCoordinate(jNode, 0))
          self.coordinate_Y.append(self.beam.ExtractCoordinate(jNode, 1))
          self.coordinate_Z.append(self.beam.ExtractCoordinate(jNode, 2))

          self.coordinate_X0.append(self.beam.ExtractCoordinate0(jNode, 0))
          self.coordinate_Y0.append(self.beam.ExtractCoordinate0(jNode, 1))
          self.coordinate_Z0.append(self.beam.ExtractCoordinate0(jNode, 2))

          self.displacement_X.append(self.beam.ExtractDisplacements(jNode, 0))
          self.displacement_Y.append(self.beam.ExtractDisplacements(jNode, 1))
          self.displacement_Z.append(self.beam.ExtractDisplacements(jNode, 2))


  def RestartLin(self):
      """ This function runs the restart and stores the results.
          Needs to be run after __SetLoads """

      self.beam.RunRestartLin(0)

      self.coordinate_X = []
      self.coordinate_Y = []
      self.coordinate_Z = []

      self.coordinate_X0 = []
      self.coordinate_Y0 = []
      self.coordinate_Z0 = []

      self.displacement_X = []
      self.displacement_Y = []
      self.displacement_Z = []

      for jNode in range(0, self.nPoint):
          self.coordinate_X.append(self.beam.ExtractCoordinate(jNode, 0))
          self.coordinate_Y.append(self.beam.ExtractCoordinate(jNode, 1))
          self.coordinate_Z.append(self.beam.ExtractCoordinate(jNode, 2))

          self.coordinate_X0.append(self.beam.ExtractCoordinate0(jNode, 0))
          self.coordinate_Y0.append(self.beam.ExtractCoordinate0(jNode, 1))
          self.coordinate_Z0.append(self.beam.ExtractCoordinate0(jNode, 2))

          self.displacement_X.append(self.beam.ExtractDisplacements(jNode, 0))
          self.displacement_Y.append(self.beam.ExtractDisplacements(jNode, 1))
          self.displacement_Z.append(self.beam.ExtractDisplacements(jNode, 2))
          
          
  def PrintDisplacements(self, iVertex):

    """ This function prints to screen the displacements on the nodes """
    print("\n--> Coord0({}): {:16.12f} {:16.12f} {:16.12f}".format(iVertex,
                                      self.coordinate_X0[iVertex],
                                      self.coordinate_Y0[iVertex],
                                      self.coordinate_Z0[iVertex]))

    print("--> Coord({}) : {:16.12f} {:16.12f} {:16.12f}".format(iVertex,
                                      self.coordinate_X[iVertex],
                                      self.coordinate_Y[iVertex],
                                      self.coordinate_Z[iVertex]))

    print("--> Displ({}) : {:16.12f} {:16.12f} {:16.12f}\n".format(iVertex,
                                      self.displacement_X[iVertex],
                                      self.displacement_Y[iVertex],
                                      self.displacement_Z[iVertex]))

  
  
  def PrintSensitivitiesAllLoads(self):

    """ This function prints the sensitivities of the objective functions for all the loads"""
    print("E', Nu' = (", self.beam.ExtractGradient_E(), self.beam.ExtractGradient_Nu(), ")")
    for iNode in range(0,self.nPoint):
       print("F'(",iNode,") = (", self.beam.ExtractLoadGradient(iNode,0), self.beam.ExtractLoadGradient(iNode,1), self.beam.ExtractLoadGradient(iNode,2), ")")
       
       
  def PrintSensitivityLoad(self, iNode):

      """ This function prints the sensitivities of the objective functions for the single load"""
      sensX = self.beam.ExtractLoadGradient(iNode,0)
      sensY = self.beam.ExtractLoadGradient(iNode,1)
      sensZ = self.beam.ExtractLoadGradient(iNode,2)

      print("F'(",iNode,") = (", sensX, sensY, sensZ, ")")

      return sensX, sensY, sensZ

  
  def PrintSensitivityE(self):

      """ This function prints the sensitivities of the objective functions for all the loads"""
      print("E' = ", self.beam.ExtractGradient_E())

      return self.beam.ExtractGradient_E()
  
  
  def PrintSensitivityPropDVs(self):

    """ This function prints the sensitivities of the objective functions wrt  Prop DVs"""
    grad_DVs =[]
    for iPDV in range(0,self.nPropDVs):
        grad_DVs.append(self.beam.ExtractPropGradient(iPDV))
        print("Prop DV '(",iPDV,") = (", grad_DVs[iPDV], ")")
    return grad_DVs

    #return self.beam.ExtractGradient_E() 
    
  def PrintSensitivityPropA(self):

      """ This function prints the sensitivities of the objective functions wrt to first Property """
      print("C_WB or A' = ", self.beam.ExtractGradient_A() )
      return self.beam.ExtractGradient_A()
  
  
  
  def TestSensitivityE(self, sensEres, sensECheck):

      """ This function prints the sensitivities of the objective functions for all the loads"""

      test_val = 100*abs((sensEres - sensECheck)/sensECheck)

      # Tolerance is set to 0.0001%
      if (test_val < 0.000001):
          print("--> Test E sensitivity: Error to reference {:10.7f}% (< 0.000001%) -> PASSED".format(test_val))
          return (0)
      else:
          print("--> Test E sensitivity: Error to reference {:10.7f}% (> 0.000001%) -> FAILED".format(test_val))
          return (1)

  def TestLoadSensitivity(self, iNode, sensX, sensY, sensZ, sensX_FD, sensY_FD, sensZ_FD):

      errorX = 100 * abs((sensX - sensX_FD) / sensX_FD)
      errorY = 100 * abs((sensY - sensY_FD) / sensY_FD)
      errorZ = 100 * abs((sensZ - sensZ_FD) / sensZ_FD)

    # Tolerance is set to 1E-8
      if errorX < 0.5 and errorY < 0.5 and errorZ < 0.5:
        print("--> Test Load sensitivity: Error to FD below 0.5% -> PASSED")
        return(0)
      else:
        print("--> Test Load sensitivity: Error to FD above 0.5% -> FAILED")
        return(1)

  def GetInitialCoordinates(self,iVertex):

    """ This function returns the initial coordinates of the structural beam model  """
    coordX = self.beam.ExtractCoordinate0(iVertex, 0)
    coordY = self.beam.ExtractCoordinate0(iVertex, 1)
    coordZ = self.beam.ExtractCoordinate0(iVertex, 2)

    return coordX, coordY, coordZ

  def ExtractDisplacements(self,iVertex):

    """ This function returns the initial coordinates of the structural beam model  """
    dispX = self.beam.ExtractDisplacements(iVertex, 0)
    dispY = self.beam.ExtractDisplacements(iVertex, 1)
    dispZ = self.beam.ExtractDisplacements(iVertex, 2)

    return dispX, dispY, dispZ

  def SetLowVerbosity(self):

      self.beam.SetLowVerbosity()

  def SetHighVerbosity(self):

      self.beam.SetHighVerbosity()



