#!/usr/bin/env python
#
# pyBeam, a Beam Solver
#
# Copyright (C) 2019 by the authors
# 
# File developers: Ruben Sanchez (SciComp, TU Kaiserslautern)
#                  Rocco Bombardieri (Carlos III University Madrid)
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
import os,os.path
import pyBeam

# ----------------------------------------------------------------------
#  Beam object
# ----------------------------------------------------------------------

class pyBeamSolver:
  """Description"""

  def __init__(self, file_dir, config_fileName):
    """ Description. """

    self.file_dir = file_dir
    self.Config_file = self.file_dir + '/' + config_fileName
    self.Config = {}
    
    print("\n---------------------------------------------------------------------------")
    print("|                                                                         |")
    print("| pyBeam, a Beam Solver - Release 0.1 (beta)                              |")
    print("|                       - https://github.com/rsanfer/pyBeam               |")
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
    self.beam = pyBeam.CBeamSolver()
    self.inputs = pyBeam.CInput(self.nPoint, self.nElem, self.nRBE2)


    # Sending to CInput object
    pyConfig.parseInput(self.Config, self.inputs, self.Constr, self.nConstr)
    # Assigning input values to the input object in C++
    self.inputs.SetParameters()
    # Initialize the input in the beam solver
    self.beam.InitializeInput(self.inputs)
    
    # Assigning values to the CNode objects in C++
    self.node = []
    for i in range(self.nPoint):
        self.node.append(pyBeam.CNode(self.node_py[i].GetID()))
        for j in range(self.nDim):
            self.node[i].InitCoordinate(j, float(self.node_py[i].GetCoord()[j][0]))
        self.beam.InitializeNode(self.node[i], i)

    # Assigning property values to the property objects in C++
    self.beam_prop = []
    for i in range(self.nProp):
        self.beam_prop.append(pyBeam.CProperty(i))
        if self.Prop[i].GetFormat() == "N":
            self.beam_prop[i].SetSectionProperties(self.Prop[i].GetA(), self.Prop[i].GetIyy(), self.Prop[i].GetIzz(), self.Prop[i].GetJt())
        elif self.Prop[i].GetFormat() == "S":
            self.beam_prop[i].SetSectionProperties(self.Prop[i].GetC_wb(), self.Prop[i].Geth(), self.Prop[i].Gett_sk(), self.Prop[i].Gett_sp(),\
            self.Prop[i].GetA_fl(), self.Prop[i].Getn_stiff(),self.Prop[i].GetA_stiff() )
            print(self.Prop[i].GetA())
        else: 
            raise ValueError("Unknown paramter for Property CARD input. Execution aborted")     
            
    # Assigning element values to the element objects in C++
    self.element = []
    for i in range(self.nElem):
        self.element.append(pyBeam.CElement(i))
        self.element[i].Initializer(self.node[self.elem_py[i].GetNodes()[0, 0] - 1], self.node[self.elem_py[i].GetNodes()[1, 0] - 1],
                               self.beam_prop[self.elem_py[i].GetProperty() - 1], self.inputs, self.elem_py[i].GetAuxVector()[0, 0],
                               self.elem_py[i].GetAuxVector()[1, 0], self.elem_py[i].GetAuxVector()[2, 0])
        self.beam.InitializeElement(self.element[i], i)

    # Here we need to pass the AeroPoint matrix of the wing grid
    # IF ANY, assigning RBE2_element values to the RBE2 objects in C++
    if self.nRBE2 != 0:
        self.RBE2 = []
        for i in range(self.nRBE2):
            self.RBE2.append(pyBeam.CRBE2(i))
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

  def GetDesignVariables(self):

      self.beam_prop = []
      for i in range(self.nProp):
          if self.Prop[i].GetFormat() == "S":
              C_wb = self.Prop[i].GetC_wb()
              h = self.Prop[i].Geth()
              A_fl = self.Prop[i].GetA_fl()
              A_stiff = self.Prop[i].GetA_stiff()
              n_stiff = self.Prop[i].Getn_stiff()
              t_sk = self.Prop[i].Gett_sk()
              t_sp = self.Prop[i].Gett_sp()


              print("C_wb =",C_wb,"h =",h,"A_fl =",A_fl,"A_stiff =",A_stiff,"n_stiff =",n_stiff,"t_sk =",t_sk,"t_sp =",t_sp)
              return C_wb,h,A_fl,A_stiff,n_stiff,t_sk,t_sp




  def CheckNewDesign(self,index_iter,C_wb_new,h_new,t_sk_new,t_sp_new,A_fl_new,n_stiff_new,A_stiff_new):
      save_path_design = '/home/marco/pyBeam/pyBeam/temp_testcase/Long_Beam_no_stiff _optimization/Design'
      save_path_prop = '/home/marco/pyBeam/pyBeam/temp_testcase/Long_Beam_no_stiff _optimization'
      self.beam_prop = []
      for i in range(self.nProp):
          if self.Prop[i].GetFormat() == "S":
              C_wb = self.Prop[i].GetC_wb()
              h = self.Prop[i].Geth()
              A_fl = self.Prop[i].GetA_fl()
              A_stiff = self.Prop[i].GetA_stiff()
              n_stiff = self.Prop[i].Getn_stiff()
              t_sk = self.Prop[i].Gett_sk()
              t_sp = self.Prop[i].Gett_sp()
      if (C_wb_new == C_wb and h_new==h and A_fl_new == A_fl and A_stiff_new == A_stiff and\
          n_stiff_new == n_stiff and t_sk_new == t_sk and t_sp_new == t_sp):
          print("Same of the previous design ")
      else:
          Name_design = os.path.join(save_path_design, str(index_iter) + "Design.txt")
          Name_prop = os.path.join(save_path_prop,"property.prt")
          with open(Name_design, "w") as f_design:
              f_design.write( str(C_wb_new) + " " + str(h_new) + " " + str(t_sk_new) + " " + str(t_sp_new) + " " + str(A_fl_new) + " " + str(
                  n_stiff_new) + " " + str(A_stiff_new))
              f_design.close()
          with open(Name_prop,"w") as f_prop:
              f_prop.write( "% ChordofWB  HeightofWB thickskin  thicksparweb Areasparcap numberstiff Astiff" +
                           " \n NPROPS=" + str(self.nProp) + "\n" + "S\n" + str(C_wb_new) + "  " + str(h_new) + "  " + str(t_sk_new) + "  " + str(
                           t_sp_new) + "  " + str(A_fl_new) + "  " + str(n_stiff_new) + "  " + str(A_stiff_new))
              f_prop.close()




  def SetLoads(self, iVertex, loadX, loadY, loadZ):

    """ This function sets the load  """
    self.beam.SetLoads(iVertex, 0, loadX)
    self.beam.SetLoads(iVertex, 1, loadY)
    self.beam.SetLoads(iVertex, 2, loadZ)

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

  def ComputeObjectiveFunction(self, iNode):

      """ This function computes the objective function (Important to be recorded) """
      displacement = self.beam.OF_NodeDisplacement(iNode)
      print("Objective Function - Displacement(", iNode, ") = ", displacement)

      return displacement

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

  def PrintSolution(self, iVertex):
    
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

    return self.coordinate_X[iVertex], self.coordinate_Y[iVertex], self.coordinate_Z[iVertex]

  def PrintPosition(self, iVertex):

      """ This function prints to screen the displacements on the nodes """
      print("--> Coord({}) : {:16.12f} {:16.12f} {:16.12f}".format(iVertex,
                                                                   self.coordinate_X[iVertex],
                                                                   self.coordinate_Y[iVertex],
                                                                   self.coordinate_Z[iVertex]))

  def PrintDisplacement(self, iVertex):

      """ This function prints to screen the displacements on the nodes """
      print("--> Displ({}) : {:16.12f} {:16.12f} {:16.12f}".format(iVertex,
                                                                   self.displacement_X[iVertex],
                                                                   self.displacement_Y[iVertex],
                                                                   self.displacement_Z[iVertex]))

  def Debug_Print(self, iElement):
      """ This function prints some input information for debugging purposes """
      self.beam.Debug_Print(iElement)


  def TestNodePosition(self, coorX, coorY, coorZ, coorX_ref, coorY_ref, coorZ_ref, tol):
    
    test_val = np.sqrt((coorX-coorX_ref)**2+
                       (coorY-coorY_ref)**2+
                       (coorZ-coorZ_ref)**2)

    if (test_val < tol):
      print("--> Test node position differs {:16.12E} from reference (< {:16.12E}) -> PASSED".format(test_val, tol))
      return(0)
    else:
      print("--> Test node position differs {:16.12E} from reference (> {:16.12E}) -> FAILED".format(test_val, tol))
      return(1)

  def SetLowVerbosity(self):

      self.beam.SetLowVerbosity()

  def SetHighVerbosity(self):

      self.beam.SetHighVerbosity()







