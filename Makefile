#============================  LIBS & INCLUDES

#EIGEN_PATH = /home/rsanfer/Software/eigen3.2.4/
EIGEN_PATH = ./pyBeam/externals/Eigen
CODI_PATH = ./pyBeam/externals/CoDiPack/include

INCLPATH =  -I$(EIGEN_PATH) -I$(CODI_PATH) $(PYTHON_INCLUDE)
#==========================================================

_SRC = ./pyBeam/src/main.cpp \
       ./pyBeam/src/ad.cpp \
       ./pyBeam/src/input.cpp \
       ./pyBeam/src/geometry.cpp \
       ./pyBeam/src/element.cpp \
       ./pyBeam/src/structure.cpp \
       ./pyBeam/src/rotations.cpp \
       ./pyBeam/src/beam.cpp

_OBJS = ./pyBeam/obj/main.o \
        ./pyBeam/obj/ad.o \
        ./pyBeam/obj/input.o \
        ./pyBeam/obj/geometry.o \
        ./pyBeam/obj/element.o \
        ./pyBeam/obj/structure.o \
        ./pyBeam/obj/rotations.o \
        ./pyBeam/obj/beam.o 

primal:
	swig -c++ -python -Wall ./pyBeam/swig/pyBeam.i
	# I need to add all the requested .cpp files
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -fPIC ./pyBeam/swig/pyBeam_wrap.cxx $(_SRC) $(INCLPATH)
	# Move objects
	mv *.o ./pyBeam/obj
	# That need to be used here so that the relative dynamic library compiles 
	g++ -O2 -shared -fPIC $(_OBJS) ./pyBeam/obj/pyBeam_wrap.o -o ./pyBeam/lib/_pyBeam.so  $(INCLPATH) 

forward:
	swig -c++ -python -Wall ./pyBeam/swig/pyBeamFM.i
	# I need to add all the requested .cpp files
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -DCODI_FORWARD_TYPE -fPIC ./pyBeam/swig/pyBeamFM_wrap.cxx $(_SRC) $(INCLPATH)
	# Move objects
	mv *.o ./pyBeam/obj	
	# That need to be used here so that the relative dynamic library compiles 
	g++ -O2 -shared -DCODI_FORWARD_TYPE -fPIC $(_OBJS) ./pyBeam/obj/pyBeamFM_wrap.o -o ./pyBeam/lib/_pyBeamFM.so  $(INCLPATH) 
	
reverse:
	swig -c++ -python -Wall ./pyBeam/swig/pyBeamAD.i
	# I need to add all the requested .cpp files
	g++ -O2 -c -w -std=gnu++11 -DLINUX=1 -DCODI_REVERSE_TYPE -fPIC ./pyBeam/swig/pyBeamAD_wrap.cxx $(_SRC) $(INCLPATH) 
	# Move objects
	mv *.o ./pyBeam/obj	
	# That need to be used here so that the relative dynamic library compiles 
	g++ -O2 -shared -DCODI_REVERSE_TYPE -fPIC $(_OBJS) ./pyBeam/obj/pyBeamAD_wrap.o -o ./pyBeam/lib/_pyBeamAD.so  $(INCLPATH) 	
	
all: primal forward reverse

RMc = rm  # remove option
# This cleans all the old object files (which are intermediate files)
clean:
	@echo "NOW CLEANING ALL" 
	-$(RMc) -f *.o
	-$(RMc) -f *.so
	-$(RMc) -f *.pyc
	-$(RMc) -f pyMLS_Cpp.py
	-$(RMc) -f *.cxx

