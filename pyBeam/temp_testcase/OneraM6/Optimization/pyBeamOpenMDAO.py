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
import math
import re
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



    def SetInitialParameters(self):

        beam = pyBeamSolverAD(self.file_dir, self.config_fileName)
        [x_tip,y_tip,z_tip] = beam.GetInitialCoordinates(self.nProp)
        [x_root, y_root, z_root] = beam.GetInitialCoordinates(0)


        #L = math.sqrt((x_tip-x_root)**2 + (y_tip-y_root)**2 + (z_tip-z_root)**2) # initial length
        L = x_tip-x_root
        print("Lunghezza ", L)

        prop_file = open('property.prt', 'r')            #read original property file
        prop_file = prop_file.read()
        dvs = re.findall(r"[-+]?\d*\.\d+|\d+", prop_file)
        dvs = np.delete(dvs, [0, 1, 2, 3, 4, 5])         # Delete old format lines from property.prt file
        DVs= np.array(list(map(float,dvs)))
        nDVs =len(DVs[0:7])                              # new format



        return DVs,nDVs,int(self.nProp),L

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


        print("lunghezza",len(DVs))
        if len(DVs) == 7:

            """ One Property"""
            Name_prop = os.path.join(final_directory, "property.prt")

            with open(Name_prop, "w") as f_prop:
                f_prop.write("% 1ChordofWB  HeightofWB thickskin  thicksparweb Areasparcap numberstiff Astiff" +
                             " \n NPROPS=" + str(self.nProp) + "\n" + "S\n" + str(DVs[0]) + "  " + str(
                    DVs[1]) + "  " + str(
                    DVs[2]) + "  " + str(DVs[3]) + "  " + str(DVs[4]) + "  " + str(int(DVs[5])) + "  " + str(DVs[6]))
                f_prop.close()

        else:

            DVs = DVs.reshape(self.nProp, 7)
            DVs = np.matrix(DVs)
            Name_prop = os.path.join(final_directory, "property_opt.prt")
            Name_prop_opt = os.path.join(final_directory, "property.prt")
            n_stiff = int(DVs[0, 5])

            with open(Name_prop, "w") as f_prop:
                for line in DVs:
                    np.savetxt(f_prop, line, fmt='%.2f')

            f_prop.close()


            with open(Name_prop, 'r') as f:
                with open(Name_prop_opt, 'w') as outfile:
                    outfile.write("% ChordofWB  HeightofWB thickskin  thicksparweb Areasparcap numberstiff Astiff" +
                            " \n NPROPS=" + str(self.nProp) + "\n" + "S\n")
                with open(Name_prop_opt, 'a') as outfile:
                    for line in f:
                        line = line.split()
                        new_line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(line[0], line[1],line[2],line[3],line[4],n_stiff,line[6])
                        outfile.write(new_line + "\n")
            os.remove(Name_prop)
            f.close()


        return final_directory

    def ComputeResponseKSStress_opt_lin(self, DVs):

        final_directory = self.NewDesign( DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        beam.SetLoads(self.Loads[0],self.Loads[1],self.Loads[2], self.Loads[3])
        beam.RunLin()
        KS=beam.ComputeResponseKSStress()

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

    def ComputeResponseKSBuckling_opt_lin(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])
        beam.RunLin()
        KS_buckl = beam.ComputeResponseKSBuckling()

        del beam.file_dir
        Name_KS_buckl = os.path.join(final_directory, "constraint_KS_buckling.txt")

        with open(Name_KS_buckl, "w") as f_KS:
            f_KS.write(str(KS_buckl))
            f_KS.close()
        print("Buckling", KS_buckl)
        return KS_buckl

    def ComputeResponseKSBuckling_opt_non_lin(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])
        beam.Run()
        KS_buckl = beam.ComputeResponseKSBuckling()
        del beam.file_dir
        Name_KS_buckl = os.path.join(final_directory, "constraint_KS_Buckling.txt")

        with open(Name_KS_buckl, "w") as f_KS:
            f_KS.write(str(KS_buckl))
            f_KS.close()

        return KS_buckl




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

    def ComputeAdjointKSstresses_opt_lin(self, DVs):

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
        beam.StopRecordingKSstresses()

        # ---------------
        beam.ComputeAdjointKSStresses()

        KSstress_sens = beam.PrintSensitivityPropDVs()

        del beam.file_dir

        return KSstress_sens




    def ComputeAdjointKSstresses_opt_non_lin(self, DVs):

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
        beam.StopRecordingKSstresses()

        # ---------------
        beam.ComputeAdjointKSStresses()
        # ----------------------

        KSstress_sens=beam.PrintSensitivityPropDVs()

        del beam.file_dir


        return KSstress_sens

    def ComputeAdjointKSBuckling_opt_lin(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        """ This function computes the response weight of the structure (important to be recorded) """

        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])
        beam.RunLin()
        # ---- Ad-hoc for

        restart_file = 'restart.pyBeam'
        restart_dir = os.getcwd()
        restart_file_location = os.path.join(restart_dir, restart_file)

        solution_file = 'solution.pyBeam'
        solution_dir = os.getcwd()
        solution_file_location = os.path.join(solution_dir, solution_file)

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

        beam.ComputeResponseKSBuckling()
        # beam.ComputeResponseSigmaBoom()

        # ---- Ad-hoc for
        # beam.StopRecordingSigmaBoom()

        # ---------------
        # beam.ComputeAdjointSigmaBoom()
        # ----------------------
        # beam.PrintSensitivitiesAllLoads()
        beam.StopRecordingKSbuckling()

        # ---------------
        beam.ComputeAdjointKSBuckling()

        KSbuckl_sens = beam.PrintSensitivityPropDVs()

        del beam.file_dir

        return KSbuckl_sens

    def ComputeAdjointKSBuckling_opt_non_lin(self, DVs):

        final_directory = self.NewDesign(DVs)
        beam = pyBeamSolverAD(final_directory, self.config_fileName)
        """ This function computes the response weight of the structure (important to be recorded) """

        beam.SetLoads(self.Loads[0], self.Loads[1], self.Loads[2], self.Loads[3])
        beam.Run()
        # ---- Ad-hoc for

        restart_file = 'restart.pyBeam'
        restart_dir = os.getcwd()
        restart_file_location = os.path.join(restart_dir, restart_file)

        solution_file = 'solution.pyBeam'
        solution_dir = os.getcwd()
        solution_file_location = os.path.join(solution_dir, solution_file)

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
        beam.ComputeResponseKSBuckling()

        # ---- Ad-hoc for
        beam.StopRecordingKSbuckling()

        # ---------------
        beam.ComputeAdjointKSBuckling()
        # ----------------------

        KSbuckl_sens = beam.PrintSensitivityPropDVs()

        del beam.file_dir

        return KSbuckl_sens





