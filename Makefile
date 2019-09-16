# To add to your .bashrc

#export PYTHONPATH=$PYTHONPATH:"/YOUR_PYBEAM_LOCATION/pyBeam/pyBeam/lib"
#export PYTHONPATH=$PYTHONPATH:"/YOUR_PYBEAM_LOCATION/pyBeam/pyBeam/pylib"
#export PYTHONPATH=$PYTHONPATH:"/YOUR_PYBEAM_LOCATION/pyBeam/pyMLS/lib"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/YOUR_PYBEAM_LOCATION/pyBeam/pyMLS/lib"
#export PYTHONPATH=$PYTHONPATH:"/YOUR_PYBEAM_LOCATION/pyBeam/pyMLS/pylib"
#export PYBEAM_INCLUDEPY="$(python-config --includes)"

#============================  LIBS & INCLUDES

# pyBeam dependencies
EIGEN_PATH = ./pyBeam/externals/Eigen
CODI_PATH = ./pyBeam/externals/CoDiPack/include

INCLPATH =  -I$(EIGEN_PATH) -I$(CODI_PATH) $(PYBEAM_INCLUDEPY)

# pyMLS dependencies
IGL_PATH = ./pyMLS/externals/libigl
ANN_PATH = ./pyMLS/externals/ann
ANN_LIB = ./pyMLS/externals/ann/lib

INCLPATH_MLS =  -I$(EIGEN_PATH) -I$(ANN_PATH) -I$(IGL_PATH) $(PYBEAM_INCLUDEPY)
LIBS_MLS = -L$(ANN_LIB) -lANN
#==========================================================

_SRC = ./pyBeam/src/main.cpp \
       ./pyBeam/src/ad.cpp \
       ./pyBeam/src/input.cpp \
       ./pyBeam/src/geometry.cpp \
       ./pyBeam/src/element.cpp \
       ./pyBeam/src/rigid_element.cpp \
       ./pyBeam/src/structure.cpp \
       ./pyBeam/src/rotations.cpp \
       ./pyBeam/src/beam.cpp \
       ./pyBeam/src/property.cpp

_SRC_MLS = ./pyMLS/src/interface.cpp

_OBJS = ./pyBeam/obj/main.o \
        ./pyBeam/obj/ad.o \
        ./pyBeam/obj/input.o \
        ./pyBeam/obj/geometry.o \
        ./pyBeam/obj/element.o \
        ./pyBeam/obj/rigid_element.o \
        ./pyBeam/obj/structure.o \
        ./pyBeam/obj/rotations.o \
        ./pyBeam/obj/beam.o \
        ./pyBeam/obj/property.o 

_OBJ_MLS = ./pyMLS/obj/interface.o \
           ./pyMLS/obj/pyMLS_wrap.o


primal:
	swig -c++ -python -Wall ./pyBeam/swig/pyBeam.i
	# Compile the objects and generate the python wrapped functions
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -fPIC ./pyBeam/swig/pyBeam_wrap.cxx $(_SRC) $(INCLPATH)
	# Move objects
	mv *.o ./pyBeam/obj
	# Compile the dynamic library
	g++ -O2 -shared -fPIC $(_OBJS) ./pyBeam/obj/pyBeam_wrap.o -o ./pyBeam/lib/_pyBeam.so  $(INCLPATH)

forward:
	swig -c++ -python -Wall ./pyBeam/swig/pyBeamFM.i
	# Compile the objects and generate the python wrapped functions
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -DCODI_FORWARD_TYPE -fPIC ./pyBeam/swig/pyBeamFM_wrap.cxx $(_SRC) $(INCLPATH)
	# Move objects
	mv *.o ./pyBeam/obj
	# Compile the dynamic library
	g++ -O2 -shared -DCODI_FORWARD_TYPE -fPIC $(_OBJS) ./pyBeam/obj/pyBeamFM_wrap.o -o ./pyBeam/lib/_pyBeamFM.so  $(INCLPATH)
	
reverse:
	swig -c++ -python -Wall ./pyBeam/swig/pyBeamAD.i
	# Compile the objects and generate the python wrapped functions
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -DCODI_REVERSE_TYPE -fPIC ./pyBeam/swig/pyBeamAD_wrap.cxx $(_SRC) $(INCLPATH)
	# Move objects
	mv *.o ./pyBeam/obj
	# Compile the dynamic library
	g++ -O2 -shared -DCODI_REVERSE_TYPE -fPIC $(_OBJS) ./pyBeam/obj/pyBeamAD_wrap.o -o ./pyBeam/lib/_pyBeamAD.so  $(INCLPATH)
	
mls_interface:
	swig -c++ -python -Wall ./pyMLS/swig/pyMLS.i
	# Compile the objects and generate the python wrapped functions
	g++ -O2 -c -fPIC -Wno-ignored-attributes -Wno-deprecated-declarations ./pyMLS/swig/pyMLS_wrap.cxx $(_SRC_MLS) $(INCLPATH_MLS)
	# Move objects
	mv *.o ./pyMLS/obj
	# Compile the dynamic library
	g++ -O2 -shared -fPIC $(_OBJ_MLS) -o ./pyMLS/lib/_pyMLS.so $(LIBS_MLS)

all: primal reverse mls_interface

RMc = rm  # remove option
# This cleans all the old object files (which are intermediate files)
clean:
	@echo "NOW CLEANING ALL" 
	-$(RMc) -f *.o
	-$(RMc) -f *.so
	-$(RMc) -f *.pyc
	-$(RMc) -f pyMLS.py
	-$(RMc) -f *.cxx

