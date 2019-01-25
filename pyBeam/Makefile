#============================  LIBS & INCLUDES

EIGEN_PATH = /home/rsanfer/Software/eigen3.2.4/
CODI_PATH = /home/rsanfer/pyBeam/pyBeam/CoDiPack/include

INCLPATH =  -I$(EIGEN_PATH) -I$(CODI_PATH)
#==========================================================


all:
	swig -c++ -python -Wall -py3 pyBeam.i
	# I need to add all the requested .cpp files used in the DLM.cpp for the relative .o files... -DDEBG
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -fPIC pyBeam_wrap.cxx ./src/ROMFM_FULLORDER.cpp ./src/input.cpp ./src/FiniteElement.cpp ./src/StructSyst.cpp ./src/Rotations.cpp ./src/beam.cpp -I/usr/include/python3.7m $(INCLPATH)
	# ... That need to be used here so that the relative pyDLM_Cpp_wrap.o compiles the correct _pyDLM_Cpp.so
	g++ -O2 -shared -fPIC ROMFM_FULLORDER.o input.o FiniteElement.o StructSyst.o Rotations.o beam.o pyBeam_wrap.o -o _pyBeam.so  $(INCLPATH) 




# (alternative version of line 2)
#g++ -O2 -c -std=gnu++11 -fPIC ./dlm/DLM.cpp ./dlm/Param.cpp $(INCLPATH) 
#g++ -O2 -c -std=gnu++11 -fPIC pyDLM_Cpp_wrap.cxx -I/usr/include/python2.7 $(INCLPATH)

#error: ImportError: './_pyDLM_Cpp.so: undefined symbol: _ZN5Param16read_AerogeneralEv'

RMc = rm  # remove option
# This cleans all the old object files (which are intermediate files)
clean:
	@echo "NOW CLEANING ALL" 
	-$(RMc) -f *.o
	-$(RMc) -f *.so
	-$(RMc) -f *.pyc
	-$(RMc) -f pyMLS_Cpp.py
	-$(RMc) -f *.cxx

