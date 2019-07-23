#!/usr/bin/env python

## \file fsi_computation.py
#  \brief Python wrapper code for FSI computation by coupling pyBeam and SU2.
#  \author David Thomas, Rocco Bombardieri, Ruben Sanchez
#  \version 7.0.0
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import time as timer

from optparse import OptionParser  # use a parser for configuration

from libFSI import FSI_config as io  # imports FSI python tools
from libFSI import Interface as FSI # imports FSI python tools
from libFSI import pyBeam_Interface as pyBeam_Interface
from libFSI import Spline_Module_class as Spline_Module

# imports the CFD (SU2) module for FSI computation
import pysu2
import pyBeam


# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():
    # --- Get the FSI config file name form the command line options --- #
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("--parallel", action="store_true",
                      help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

    (options, args) = parser.parse_args()

    if options.with_MPI:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        numberPart = comm.Get_size()
        have_MPI = True
    else:
        comm = 0
        myid = 0
        numberPart = 1
        have_MPI = False

    rootProcess = 0

    # --- Set the working directory --- #
    if myid == rootProcess:
        if os.getcwd() not in sys.path:
            sys.path.append(os.getcwd())
            print("Setting working directory : {}".format(os.getcwd()))
        else:
            print("Working directory is set to {}".format(os.getcwd()))

    # starts timer
    start = timer.time()

    confFile = str(options.filename)

    FSI_config = io.FSIConfig(confFile)  # FSI configuration file
    CFD_ConFile = FSI_config['SU2_CONFIG']  # CFD configuration file
    CSD_ConFile = FSI_config['PYBEAM_CONFIG']  # CSD configuration file
    MLS_confFile = FSI_config['MLS_CONFIG_FILE_NAME']  # MLS configuration file

    if have_MPI:
        comm.barrier()

    # --- Initialize the fluid solver: SU2 --- #
    if myid == rootProcess:
        print('\n***************************** Initializing SU2 **************************************')
    try:
        FluidSolver = pysu2.CSinglezoneDriver(CFD_ConFile, 1, FSI_config['NDIM'], comm)
    except TypeError as exception:
        print('A TypeError occured in pysu2.CSingleZoneDriver : ', exception)
        if have_MPI:
            print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
        else:
            print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
        return

    if have_MPI:
        comm.barrier()

    # --- Initialize the solid solver: pyBeam --- #
    if myid == rootProcess:
        print('\n***************************** Initializing pyBeam ************************************')
        try:
            SolidSolver = pyBeam_Interface.pyBeamSolver(CSD_ConFile)
        except TypeError as exception:
            print('ERROR building the Solid Solver: ', exception)
    else:
        SolidSolver = None

    if have_MPI:
        comm.barrier()

    # --- Initialize and set the coupling environment --- #
    if myid == rootProcess:
        print('\n***************************** Initializing FSI interface *****************************')
    try:
        FSIInterface = FSI.Interface(FSI_config, FluidSolver, SolidSolver, None, have_MPI)
    except TypeError as exception:
        print('ERROR building the FSI Interface: ', exception)

    if have_MPI:
        comm.barrier()


    if myid == rootProcess:
        print('\n***************************** Connect fluid and solid solvers *****************************')
    try:
        FSIInterface.connect(FSI_config, FluidSolver, SolidSolver)
    except TypeError as exception:
        print('ERROR building the Interpolation Interface: ', exception)

    if have_MPI:
        comm.barrier()

    if myid == rootProcess:  # we perform this calculation on the root core
        print('\n***************************** Initializing MLS Interpolation *************************')
        try:
            MLS = Spline_Module.MLS_Spline(MLS_confFile, FSIInterface.nDim,
                                          FSIInterface.globalFluidCoordinates, FSIInterface.globalSolidCoordinates,
                                          FSI_config)
        except TypeError as exception:
            print('ERROR building the MLS Interpolation: ', exception)

    else:
        MLS = None

    if have_MPI:
        comm.barrier()

    # Run the solver
    if myid == 0:
        print("\n------------------------------ Begin Solver -----------------------------\n")
    sys.stdout.flush()
    if options.with_MPI:
        comm.Barrier()

    FSIInterface.SteadyFSI(FSI_config, FluidSolver, SolidSolver, MLS)

    # Postprocess the solver and exit cleanly
    FluidSolver.Postprocessing()

    if FluidSolver is not None:
        del FluidSolver

    #
    # # --- Initialize and set the coupling environment --- #
    # if myid == rootProcess:
    #     print('\n***************************** Initializing FSI interface *****************************')
    # try:
    #     FSIInterface = FSI.Interface(FSI_config, FluidSolver, SolidSolver, MLS, have_MPI)
    # except TypeError as exception:
    #     print('ERROR building the FSI Interface: ', exception)
    #
    # if have_MPI == True:
    #     comm.barrier()
    #
    #
    # if myid == rootProcess:
    #     print('\n***************************** Mapping fluid-solid interfaces *****************************')
    # try:
    #     FSIInterface.interfaceMapping(FluidSolver, SolidSolver, FSI_config)
    # except TypeError as exception:
    #     print('ERROR building the Interface Mapping: ', exception)
    #
    # # --- Launch a steady or unsteady FSI computation --- #
    # try:
    #     FSIInterface.SteadyFSI(FSI_config, FluidSolver, SolidSolver, MLS)
    # except NameError as exception:
    #     if myid == rootProcess:
    #         print('An NameError occured in FSIInterface.SteadyFSI : ', exception)
    # except TypeError as exception:
    #     if myid == rootProcess:
    #         print('A TypeError occured in FSIInterface.SteadyFSI : ', exception)
    # except KeyboardInterrupt as exception:
    #     if myid == rootProcess:
    #         print('A KeyboardInterrupt occured in FSIInterface.SteadyFSI : ', exception)
    #
    # if have_MPI == True:
    #     comm.barrier()
    #
    # # --- Exit cleanly the fluid and solid solvers --- #
    # FluidSolver.Postprocessing()
    # if myid == rootProcess:
    #     SolidSolver.exit()
    #
    # if have_MPI == True:
    #     comm.barrier()
    #
    # # stops timer
    # stop = timer.time()
    # elapsedTime = stop - start
    #
    # if myid == rootProcess:
    #     print("\n Computation successfully performed in {} seconds.".format(elapsedTime))

    return


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
