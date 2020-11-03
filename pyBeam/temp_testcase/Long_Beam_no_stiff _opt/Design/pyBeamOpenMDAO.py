#!/usr/bin/env python3
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

from pyBeamIO import pyBeamConfig as pyConfig
from pyBeamIO import pyBeamInput as pyInput
import numpy as np
import os,os.path
import shutil
import pyBeam
from pyBeamLibAD import pyBeamSolverAD



class pyBeamOpt:
    def __init__(self, file_dir, config_fileName,Loads_Values):
        self.file_dir = file_dir
        self.config_fileName = config_fileName
        self.Config_file = self.file_dir + '/' + config_fileName
        self.Config = {}
        self.Config = pyConfig.pyBeamConfig(self.Config_file)  # Beam configuration file
        self.Property = self.file_dir + '/' + self.Config['PROPERTY_FILE']
        # Parsing Property file
        self.Prop, self.nProp = pyInput.readProp(self.Property)
        self.Loads = Loads_Values

    def NewDesign(self,DVs):
        """ This function create a folder in which are copied the input files,appending the new Design Variables (DVs) in property file   """
        path = os.path.abspath(os.path.join(self.file_dir, ".."))
        final_directory = os.path.join(path, r'Optimization')
        if not os.path.exists(final_directory):
            os.makedirs(final_directory)

        file_names = os.listdir(self.file_dir)

        for file_name in file_names:
            new_path = os.path.join(self.file_dir, file_name)
            if os.path.isdir(new_path):
                continue
            shutil.copy(os.path.join(self.file_dir, file_name), final_directory)

        Name_prop = os.path.join(final_directory, "property.prt")

        with open(Name_prop, "w") as f_prop:
            f_prop.write("% ChordofWB  HeightofWB thickskin  thicksparweb Areasparcap numberstiff Astiff" +
                         " \n NPROPS=" + str(self.nProp) + "\n" + "S\n" + str(DVs[0]) + "  " + str(
                DVs[1]) + "  " + str(
                DVs[2]) + "  " + str(DVs[3]) + "  " + str(DVs[4]) + "  " + str(DVs[5]) + "  " + str(DVs[6]))
            f_prop.close()

        return final_directory

    def ComputeResponseKSStress_opt_lin(self, DVs):

        final_directory = self.NewDesign( DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        beam.SetLoads(self.Loads[0],self.Loads[1],self.Loads[2], self.Loads[3])
        beam.RunLin()
        KS=beam.ComputeResponseKSStress()
        #beam.ComputeResponseSigmaBoom()

        del beam.file_dir
        Name_KS = os.path.join(final_directory, "constraint_KS.txt")

        with open(Name_KS, "w") as f_KS:
            f_KS.write(str(KS))
            f_KS.close()

        return KS

    def ComputeResponseKSStress_opt_non_lin(self, DVs):

        final_directory = self.NewDesign( DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        beam.SetLoads(self.Loads[0],self.Loads[1],self.Loads[2], self.Loads[3])
        beam.Run()
        KS= beam.ComputeResponseKSStress()
        del beam.file_dir
        Name_KS = os.path.join(final_directory, "constraint_KS.txt")

        with open(Name_KS, "w") as f_KS:
            f_KS.write(str(KS))
            f_KS.close()

        return KS




    def ComputeWeight_opt( self,  DVs):

          final_directory = self.NewDesign(DVs)
          beam = pyBeamSolverAD(final_directory, self.config_fileName)
          """ This function computes the response weight of the structure (important to be recorded) """
          weight = beam.ComputeWeight()
          del beam.file_dir

          Name_weight = os.path.join(final_directory,"weight.txt")

          with open(Name_weight, "w") as f_weight:
              f_weight.write(str(weight))
              f_weight.close()

          return weight

    def ComputeAdjointWeight_opt(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        """ This function computes the response weight of the structure (important to be recorded) """

        # ---- Ad-hoc for
        beam.StartRecording()
        beam.SetDependencies()
        # ----------------------
        beam.ComputeWeight()

        beam.StopRecordingWeight()

        beam.ComputeAdjointWeight()

        Weight_sens = beam.PrintSensitivityPropDVs()


        del beam.file_dir


        Name_weight_sens = os.path.join(final_directory, "weight_sens.txt")

        with open(Name_weight_sens, "w") as f_weight:
            f_weight.write(str(Weight_sens))
            f_weight.close()

        return Weight_sens

    def ComputeAdjointKS_opt_lin(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        """ This function computes the response weight of the structure (important to be recorded) """

        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])
        beam.RunLin()
        # ---- Ad-hoc for

        restart_file = 'restart.pyBeam'
        restart_dir = os.getcwd()
        restart_file_location = os.path.join( restart_dir, restart_file)

        solution_file='solution.pyBeam'
        solution_dir = os.getcwd()
        solution_file_location=os.path.join(solution_dir, solution_file)

        with open(restart_file_location, 'r') as f1:
            with open(solution_file_location, 'w') as f2:
                for line in f1:
                    f2.write(line)

        beam.ReadRestart()

        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])


        beam.StartRecording()
        beam.SetDependencies()
        # ----------------------
        beam.RestartLin()

        beam.ComputeResponseKSStress()
        #beam.ComputeResponseSigmaBoom()

        # ---- Ad-hoc for
        #beam.StopRecordingSigmaBoom()

        # ---------------
        #beam.ComputeAdjointSigmaBoom()
        # ----------------------
        # beam.PrintSensitivitiesAllLoads()
        beam.StopRecordingKS()

        # ---------------
        beam.ComputeAdjointKS()

        KS_sens = beam.PrintSensitivityPropDVs()

        del beam.file_dir


        Name_KS_sens = os.path.join(final_directory, "KS_sens.txt")

        with open(Name_KS_sens, "w") as f_KS_sens:
            f_KS_sens.write(str(KS_sens))
            f_KS_sens.close()

        return KS_sens




    def ComputeAdjointKS_opt_non_lin(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        """ This function computes the response weight of the structure (important to be recorded) """

        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])
        beam.Run()
        # ---- Ad-hoc for

        restart_file = 'restart.pyBeam'
        restart_dir = os.getcwd()
        restart_file_location = os.path.join( restart_dir, restart_file)

        solution_file='solution.pyBeam'
        solution_dir = os.getcwd()
        solution_file_location=os.path.join(solution_dir, solution_file)

        with open(restart_file_location, 'r') as f1:
            with open(solution_file_location, 'w') as f2:
                for line in f1:
                    f2.write(line)

        beam.ReadRestart()

        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])


        beam.StartRecording()
        beam.SetDependencies()
        # ----------------------
        beam.Restart()
        beam.ComputeResponseKSStress()


        # ---- Ad-hoc for
        beam.StopRecordingKS()

        # ---------------
        beam.ComputeAdjointKS()
        # ----------------------

        KS_sens=beam.PrintSensitivityPropDVs()

        del beam.file_dir


        Name_KS_sens = os.path.join(final_directory, "KS_sens.txt")

        with open(Name_KS_sens, "w") as f_KS_sens:
            f_KS_sens.write(str(KS_sens))
            f_KS_sens.close()

        return KS_sens





