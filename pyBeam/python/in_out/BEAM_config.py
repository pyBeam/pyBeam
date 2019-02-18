import os, sys, shutil, copy
from util import switch
#import ast
from math import *

# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class BEAMConfig:
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
		if case("N_DOF")		      : pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK 
                if case("FOLLOWER_FLAG")		      : pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK       
		if case("LOAD_STEPS")	      : pass  # in case: CSD_SOLVER = NITRO_FRAMEWORK      
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
                if case("B_LENGTH")                 : pass
		if case("LOAD")		      : pass
                if case("CONV_CRITERIUM")          : 
		    self._ConfigContent[this_param] = float(this_value)
		    break

	        #string values  MEMO_GEN_FORCE_OUTPUT
                if case("B_MESH")          :       pass           
                if case("B_PROPERTY")                 :               
		    self._ConfigContent[this_param] = this_value
		    break
                              
 	        if case():
		    print(this_param + " is an invalid option !")
		    break