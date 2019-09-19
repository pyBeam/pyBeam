#!/usr/bin/env python
#
# pyMLS, an open-source Moving Least Squares library
#
# Copyright (C) 2019 by the authors
# 
# File developers: Rocco Bombardieri (Carlos III University Madrid)
#                  Rauno Cavallaro (Carlos III University Madrid)
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

from util.switch import switch

# ----------------------------------------------------------------------
#  MLS Configuration Class
# ----------------------------------------------------------------------

class pyMLSConfig:
    """
    Class that contains all the parameters coming from the MLS configuration file.
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
                if case("MODES"): pass
                if case("POLY"): pass
                if case("WEIGHT"): pass
                if case("POINTS"): 
                    self._ConfigContent[this_param] = int(this_value)
                    break

            #float values
                if case("MAGNIF_FACTOR"): pass
                if case("RMAX"): pass                
                if case("DELTA"): pass
                if case("TOLL_SVD"): pass
                if case("FREQ"):
                    self._ConfigContent[this_param] = float(this_value)
                    break

            #string values
                if case("STRUCTURAL_MODES_FILE_NAME"): pass
                if case("FORMAT_MODES"): pass
                if case("DEBUG"):
                    self._ConfigContent[this_param] = this_value
                    break

                if case():
                    print(this_param + " is an invalid option !")
                    break
