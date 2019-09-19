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
import os, sys, shutil, copy
from util.switch import switch
#import ast
from math import *

# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class pyBeamConfig:
    """
    Class that contains all the parameters coming from the BEAM configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self,FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.readConfig()

    def __str__(self):
        tempString = str()
        for key, value in self._ConfigContent.items():
          tempString += "{} = {}\n".format(key,value)
        return tempString

    def __getitem__(self,key):
        return self._ConfigContent[key]

    def __setitem__(self, key, value):
        self._ConfigContent[key] = value

    def readConfig(self):
        input_file = open(self.ConfigFileName)
        while 1:
            line = input_file.readline()
            if not line:
                break
            # remove line returns
            line = line.strip('\r\n')
            # make sure it has useful data
            if (not "=" in line) or (line[0] == '%'):
                continue
            # split across equal sign
            line = line.split("=",1)
            this_param = line[0].strip()
            this_value = line[1].strip()

            for case in switch(this_param):
                #integer values
                if case("RESTART")     :       pass
                if case("LOAD_STEPS")	      : pass
                if case("N_STRUCT_ITER")	                      : 
                    self._ConfigContent[this_param] = int(this_value)
                    break
                #float values
                if case("W_THICKNESS")         : pass                 
                if case("W_HEIGHT")        : pass                
                if case("F_WIDTH")                      : pass
                if case("Y_MODULUS")                      : pass
                if case("POISSON")                      : pass
                if case("RHO")                      : pass
                if case("LOAD")		      : pass
                if case("CONV_CRITERIUM")          : 
                    self._ConfigContent[this_param] = float(this_value)
                    break
                #string values  MEMO_GEN_FORCE_OUTPUT
                if case("TOLERANCE_LINSOL"):
                    self._ConfigContent[this_param] = float(this_value)
                    break
                if case("KIND_LINSOL"):
                    #print(this_value)
                    if this_value == "PartialPivLu":
                        self._ConfigContent[this_param] = 1
                        break
                    elif this_value == "FullPivLu":
                        self._ConfigContent[this_param] = 2
                        break
                    elif this_value == "HouseholderQr":
                        self._ConfigContent[this_param] = 3
                        break
                    elif this_value == "ColPivHouseholderQr":
                        self._ConfigContent[this_param] = 4
                        break
                    elif this_value == "FullPivHouseholderQr":
                        self._ConfigContent[this_param] = 5
                        break
                    elif this_value == "LLT":
                        self._ConfigContent[this_param] = 6
                        break
                    elif this_value == "LDLT":
                        self._ConfigContent[this_param] = 7
                        break
                    else:
                        self._ConfigContent[this_param] = 2
                        break
                if case("B_MESH")          :       pass
                if case("B_PROPERTY")                 :               
                    self._ConfigContent[this_param] = this_value
                    break
                if case():
                    print(this_param + " is an invalid option !")
                    break

def parseInput(BEAM_config, inputs,Constr, nConstr):
    
    inputs.SetYoungModulus(BEAM_config['Y_MODULUS'])
    inputs.SetPoisson(BEAM_config['POISSON'])
    inputs.SetDensity(BEAM_config['RHO'])
    inputs.SetLoadSteps(BEAM_config['LOAD_STEPS'])
    inputs.SetNStructIter(BEAM_config['N_STRUCT_ITER'])
    inputs.SetConvCriterium(BEAM_config['CONV_CRITERIUM'])
    inputs.SetTolerance_LinSol(BEAM_config['TOLERANCE_LINSOL'])
    inputs.SetKind_LinSol(BEAM_config['KIND_LINSOL'])
    
    #Now setting the constraints
    inputs.SetnConstr(nConstr)
    for iConstr in range(nConstr): 
        inputs.SetSingleConstr( iConstr, int(Constr[iConstr,0]), int(Constr[iConstr,1]) )