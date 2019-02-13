#============================  LIBS & INCLUDES

#EIGEN_PATH = /home/rsanfer/Software/eigen3.2.4/
EIGEN_PATH = ./externals/Eigen
CODI_PATH = ./externals/CoDiPack/include

INCLPATH =  -I$(EIGEN_PATH) -I$(CODI_PATH)
#==========================================================

_SRC = ./src/ROMFM_FULLORDER.cpp \
       ./src/ad.cpp \
       ./src/input.cpp \
       ./src/geometry.cpp \
       ./src/FiniteElement.cpp \
       ./src/StructSyst.cpp \
       ./src/Rotations.cpp \
       ./src/beam.cpp

_OBJS = ./obj/ROMFM_FULLORDER.o \
        ./obj/ad.o \
        ./obj/input.o \
        ./obj/geometry.o \
        ./obj/FiniteElement.o \
        ./obj/StructSyst.o \
        ./obj/Rotations.o \
        ./obj/beam.o 

primal:
	swig -c++ -python -Wall -py3 ./swig/pyBeam.i
	# I need to add all the requested .cpp files
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -fPIC ./swig/pyBeam_wrap.cxx $(_SRC) -I/usr/include/python3.7m $(INCLPATH)
	# Move objects
	mv *.o ./obj
	# That need to be used here so that the relative dynamic library compiles 
	g++ -O2 -shared -fPIC $(_OBJS) ./obj/pyBeam_wrap.o -o ./lib/_pyBeam.so  $(INCLPATH) 

forward:
	swig -c++ -python -Wall -py3 ./swig/pyBeamFM.i
	# I need to add all the requested .cpp files
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -DCODI_FORWARD_TYPE -fPIC ./swig/pyBeamFM_wrap.cxx $(_SRC) -I/usr/include/python3.7m $(INCLPATH)
	# Move objects
	mv *.o ./obj	
	# That need to be used here so that the relative dynamic library compiles 
	g++ -O2 -shared -DCODI_FORWARD_TYPE -fPIC $(_OBJS) ./obj/pyBeamFM_wrap.o -o ./lib/_pyBeamFM.so  $(INCLPATH) 
	
reverse:
	swig -c++ -python -Wall -py3 ./swig/pyBeamAD.i
	# I need to add all the requested .cpp files
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -DCODI_REVERSE_TYPE -fPIC ./swig/pyBeamAD_wrap.cxx $(_SRC) -I/usr/include/python3.7m $(INCLPATH)
	# Move objects
	mv *.o ./obj	
	# That need to be used here so that the relative dynamic library compiles 
	g++ -O2 -shared -DCODI_REVERSE_TYPE -fPIC $(_OBJS) ./obj/pyBeamAD_wrap.o -o ./lib/_pyBeamAD.so  $(INCLPATH) 	
	
all: primal forward reverse


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

