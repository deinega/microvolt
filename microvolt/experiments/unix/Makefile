#----------------------------------------------------------------------
#   Microvolt
#----------------------------------------------------------------------

#_OPENBLAS=true
#_PARDISO=true

#SHELL = /bin/sh

# C++ compiler (default value)
CPP = g++  

# Optimization flag (default value)
OPT_FLAG = -O2

# We use UNIX preprocessor flag for code which can be executed on UNIX only (not Windows)
PREPROC_FLAG := -DUNIX $(PREPROC_FLAG)

ifeq (true,$(_OPENBLAS))
PREPROC_FLAG := -DUSE_LAPACK $(PREPROC_FLAG)
USRLIB := $(USRLIB) -L/usr/local/lib -lopenblas
endif

ifeq (true,$(_PARDISO))
PREPROC_FLAG := -DUSE_PARDISO $(PREPROC_FLAG)
USRLIB := $(USRLIB) -L/usr/local/lib -lpardiso500-GNU472-X86-64
endif

ifeq (true,$(_MKL))
PREPROC_FLAG := -DUSE_MKL $(PREPROC_FLAG)
USRLIB := $(USRLIB) -L$(MKLROOT) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
endif

ROOTDIR  := ../..
UTILDIR  = $(ROOTDIR)/../ivutils

IVINCDIR   = $(UTILDIR)/include
INCDIR := $(ROOTDIR)/include
SRCDIR := $(ROOTDIR)/src
TSTDIR   = $(ROOTDIR)/experiments

INCLUDE_FLAG := -I$(IVINCDIR) -I$(INCDIR)

# VPATH is list of directories where make looks for source files.
VPATH	= $(TSTDIR):$(UTILDIR)/src:$(SRCDIR)

SRC_UTILS = region_2.cpp region_3.cpp gnudump.cpp grid.cpp vector_set.cpp \
	logexc.cpp component.cpp string_utils.cpp ifmpi.cpp \
	grid_t.cpp contour_t.cpp linsolv.cpp

SRC_TRANS = sc_mediumbase.cpp sc_mediumfactory.cpp sc_alloyfactory.cpp \
	sc_generation.cpp sc_dump.cpp sc_fd.cpp sc_mesh.cpp sc_fix.cpp \
	sc_detector.cpp sc_ui.cpp 

SRC_MAIN = diode.cpp 
#SRC_MAIN = multijunction.cpp 
#SRC_MAIN = nanowire.cpp 
#SRC_MAIN = texture.cpp 

# list of compiled source .cpp files
SRC = $(SRC_UTILS) $(SRC_TRANS) $(SRC_MAIN)

# list of resulted object .o files
OBJ = $(SRC:.cpp=.o)
	
# Make the default target - executable file 
.PHONY: all
all:: transport.x;

# suffix rules
.SUFFIXES : .o .cpp 

# rules to get .o targets from .cpp prerequisites.
# $< - first prerequisite of the rule, used to put source .cpp file after -c
.cpp.o :
	$(CPP) $(OPT_FLAG) -MD $(PREPROC_FLAG) $(INCLUDE_FLAG) -c $<

# Target rules for executable specify prerequisite objects and linking command.
transport.x: $(OBJ)
	$(CPP) $(OPT_FLAG) $(OBJ) \
	$(USRLIB) \
	-fopenmp -o microvolt.x;

# Include individual dependencies.
# Substitution reference is used, $(VAR:a=b) is VAR where all a is changed on b
# -include is include and ingore unexisting files
DEPENDS = $(OBJ:.o=.d)
-include $(DEPENDS)

clean:
	rm -rf *.d *.o microvolt.x

