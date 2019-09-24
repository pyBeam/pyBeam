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

from pyMLSIO import pyMLSConfig as io
import numpy as np
import pyMLS as Spline

# ----------------------------------------------------------------------
#  MLS_Spline Interface Class
# ----------------------------------------------------------------------

class pyMLSInterface:

    """
    MLS_Spline class that handles fluid/solid solver synchronisation and communication
    """

    def __init__(self, MLS_Config_File, AeroNodes, StructNodes):
        """
         Class constructor. Declare some variables and do some screen outputs.
        """
        # AeroPoint is numpy matrix nAeropoint*3 (AeroPoint[markers[FSI_marker])
        # BoundElem is an object. here we need the function .GetNodes() or something similar
        # The MLS configurations parameters are stored from the MLS input file
        print("Storing MLS parameters from input file ")
        MLS_conf = io.pyMLSConfig(MLS_Config_File)

        self.nAeroNodes = np.shape(AeroNodes)[0]
        self.nStructNodes = np.shape(StructNodes)[0]

        lenAeroNodes = self.nAeroNodes * 3
        lenStructNodes = self.nStructNodes * 3

        # Performing the meshless method
        print("Performing the Meshless Method")
        # Arrange structural nodes in the wrapped standard vector

        str_data_std = Spline.DoubleVector(lenStructNodes)
        l = 0
        for i in range(0, 3):
            for j in range(0, self.nStructNodes):
                str_data_std[l] = float(StructNodes[j][i])  # str_data[j,i]
                l = l + 1

        # Arrange aerodynamic nodes in the wrapped standard vector
        aero_data_std = Spline.DoubleVector(lenAeroNodes)
        l = 0
        for i in range(0, 3):
            for j in range(0, self.nAeroNodes):
                aero_data_std[l] = float(AeroNodes[j][i])
                l = l + 1

        interpolation_matrix_std = Spline.DoubleVector(self.nAeroNodes * self.nStructNodes)
        norm_err_std = Spline.DoubleVector(self.nAeroNodes)

        Spline.mls_interface(interpolation_matrix_std, norm_err_std, self.nStructNodes, self.nAeroNodes, str_data_std,
                             aero_data_std, MLS_conf['POLY'], MLS_conf['WEIGHT'], MLS_conf['POINTS'],
                             MLS_conf['RMAX'], MLS_conf['DELTA'], MLS_conf['TOLL_SVD'])

        # --- OUTPUT ----------------------------------------------------------------
        self.interpolation_matrix = np.zeros((self.nAeroNodes, self.nStructNodes))
        l = 0
        for i in range(0, self.nStructNodes):
            for j in range(0, self.nAeroNodes):
                self.interpolation_matrix[j][i] = interpolation_matrix_std[l]
                l = l + 1

        print("Splining: norm of interpolation error over nodes position = {}".format(np.linalg.norm(norm_err_std)))
