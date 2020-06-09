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

# Imports from pyBeam (need to be reviewed)
import sys, os
from pyBeamLib import pyBeamSolver as pyBeamSolver_pyBeam
from pyBeamLibAD import pyBeamSolverAD as pyBeamSolverAD_pyBeam

# ----------------------------------------------------------------------
#  Beam object
# ----------------------------------------------------------------------

class pyBeamSolver:
  """Class to handle the primal beam solver"""

  def __init__(self, config_fileName):
    """ Class initialization. """

    self.Config_file = config_fileName
    self.Config = {}

    print("\n----------- Configuring the primal solver in pyBeam ----------------")
    
    # Load testcase directory directory
    file_dir = os.getcwd()

    self.beam = pyBeamSolver_pyBeam(file_dir, self.Config_file)
    
    self.nPoint = self.beam.nPoint


  def SetLoads(self, iVertex, loadX, loadY, loadZ):

    """ This function sets the load  """
    self.beam.SetLoads(iVertex, loadX, loadY, loadZ)


  def GetInitialCoordinates(self,iVertex):

    """ This function returns the initial coordinates of the structural beam model  """
    coordX, coordY, coordZ = self.beam.GetInitialCoordinates(iVertex)

    return coordX, coordY, coordZ

  def ExtractDisplacements(self,iVertex):

    """ This function returns the initial coordinates of the structural beam model  """
    dispX, dispY, dispZ = self.beam.ExtractDisplacements(iVertex)

    return dispX, dispY, dispZ


  def run(self):
    
    """This function runs the solver. Needs to be run after SetLoads"""
 
    success =True  
    try:
       self.beam.Run()
    except:
       print("pyBeam error: check history") 
       success = False
       
    return success   
       
class pyBeamADSolver:
  """Description"""

  def __init__(self, config_fileName):
    """ Description. """

    self.Config_file = config_fileName
    self.Config = {}

    print("\n----------- Configuring the AD solver in pyBeam --------------------")

    # Load testcase directory directory
    file_dir = os.getcwd()

    self.beam = pyBeamSolverAD_pyBeam(file_dir, self.Config_file)

    # Some intermediate variables
    self.nPoint = self.beam.nPoint

    # Store the history of the sensitivity to check convergence
    self.sens_file = open("sensitivity_E.dat", "w")
    self.sens_file.write("Sensitivity E\n")
    self.sens_file.close()


  def SetLoads(self, iVertex, loadX, loadY, loadZ):

    """ This function sets the load  """
    self.beam.SetLoads(iVertex, loadX, loadY, loadZ)

  def SetDisplacementAdjoint(self, iVertex, adjX, adjY, adjZ):

    """ This function sets the adjoint cross-dependency of the displacement """
    self.beam.SetDisplacementAdjoint(iVertex, adjX, adjY, adjZ)

  def GetLoadAdjoint(self, iVertex):

    """ This function gets the adjoint of the loads """
    sensX, sensY, sensZ = self.beam.GetLoadSensitivity(iVertex)
    return sensX, sensY, sensZ

  def RecordSolver(self):
    """This function completes the pyBeam recording. Needs to be run after_SetLoads"""
    self.beam.ReadRestart()
    self.beam.StartRecording()
    self.beam.SetDependencies()
    self.beam.Restart()
    self.beam.StopRecording()

  def GetInitialCoordinates(self,iVertex):

    """ This function returns the initial coordinates of the structural beam model  """
    coordX, coordY, coordZ = self.beam.GetInitialCoordinates(iVertex)

    return coordX, coordY, coordZ

  def ExtractDisplacements(self,iVertex):

    """ This function returns the initial coordinates of the structural beam model  """
    dispX, dispY, dispZ = self.beam.ExtractDisplacements(iVertex)

    return dispX, dispY, dispZ


  def RunAdjoint(self):
    
    """This function runs the adjoint solver. Needs to be run after RecordSolver"""
    self.beam.ComputeAdjoint()
    sens_E = self.beam.PrintSensitivityE()
    self.sens_file = open("sensitivity_E.dat", "a")
    self.sens_file.write(str(sens_E) + "\n")
    self.sens_file.close()


