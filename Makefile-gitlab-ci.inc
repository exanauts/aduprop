CXX := mpicxx
CFLAGS := -fopenmp -DADUPROP_MPI -march=native -g -O2 -Wall 
# DEBUG
MATH_LIBS := -fopenmp -llapack -lblas
CODI_DIR := /builds/aduprop/aduprop/CoDiPack
ADUPROP_DIR = /builds/aduprop/aduprop 
#EIGEN_DIR := $(ADUPROP_DIR)/eigen/install/include/eigen3
echo "ADUPROP_DIR: $ADUPROP_DIR"
echo "CODI_DIR: $CODI_DIR"

# DON'T EDIT BELOW THIS LINE

ifndef CODI_DIR
 $(error Environment variable CODI_DIR is undefined)
endif

INCLUDE_DIR := /builds/aduprop/aduprop/include
SRC_DIR := /builds/aduprop/aduprop/src
echo "SRC_DIR: $SRC_DIR"
echo "INCLUDE_DIR: $INCLUDE_DIR"
CFLAGS += -std=c++11 -I$(CODI_DIR)/include -I$(INCLUDE_DIR)
#CFLAGS += -std=c++11 -I$(CODI_DIR)/include -I$(EIGEN_DIR) -I$(INCLUDE_DIR)
CFLAGS += -I/builds/aduprop/aduprop/CONTRIB
LDLIBS := -L/builds/aduprop/aduprop -laduprop $(MATH_LIBS)
AR := ar
AROPT := rcs

HEADERFILES := ad.hpp alg.hpp linsolve.hpp tensor.hpp

HEADERS := $(addprefix /builds/aduprop/aduprop/include/, $(HEADERFILES))


# HDF5 support
#ifneq ($(HDF_INSTALL),)
  CFLAGS += -I/usr/include/hdf5/serial
  CFLAGS += -DHDF5
  LDLIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
  #LDLIBS += /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a
  #LDLIBS += /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a
  LDLIBS += -lhdf5_hl -lhdf5 -lpthread -ldl -lsz -lz -lm
  #endif
