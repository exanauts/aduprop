# Compile options
# Enable MPI
ADUPROP_MPI=1
#
# Enable Eigen dense linear solver
# ADUPROP_EIGEN=1
#
# Enable Eigen sparse linear solver
# ADUPROP_EIGEN_SPARSE=1

CXX = mpicxx
CFLAGS = -fopenmp -march=native -O3 -Wall 
# DEBUG
MATH_LIBS = -fopenmp -llapack -lblas
CODI_DIR = /builds/aduprop/aduprop/CoDiPack
ADUPROP_DIR = /builds/aduprop/aduprop 
#EIGEN_DIR = $(ADUPROP_DIR)/eigen/install/include/eigen3

# DON'T EDIT BELOW THIS LINE

ifndef CODI_DIR
 $(error Environment variable CODI_DIR is undefined)
endif

ifdef ADUPROP_EIGEN
ifndef EIGEN_DIR
 $(error Environment variable EIGEN_DIR is undefined)
endif
endif
ifdef ADUPROP_EIGEN_SPARSE
ifndef EIGEN_DIR
 $(error Environment variable EIGEN_DIR is undefined)
endif
endif

INCLUDE_DIR := /builds/aduprop/aduprop/include
SRC_DIR := /builds/aduprop/aduprop/src
CFLAGS += -std=c++11 -I$(CODI_DIR)/include -I$(INCLUDE_DIR)
ifdef ADUPROP_MPI
 CFLAGS += -I$(EIGEN_DIR) -DADUPROP_MPI
endif

ifdef ADUPROP_EIGEN
 CFLAGS += -I$(EIGEN_DIR) -DADUPROP_EIGEN
endif
ifdef ADUPROP_EIGEN_SPARSE
 CFLAGS += -I$(EIGEN_DIR) -DADUPROP_EIGEN_SPARSE
endif
CFLAGS += -I/builds/aduprop/aduprop/CONTRIB
LDLIBS = -L/builds/aduprop/aduprop -laduprop $(MATH_LIBS)
AR = ar
AROPT = rcs

HEADERFILES = ad.hpp alg.hpp linsolve.hpp tensor.hpp

HEADERS = $(addprefix /builds/aduprop/aduprop/include/, $(HEADERFILES))


# HDF5 support
#ifneq ($(HDF_INSTALL),)
  CFLAGS += -I/usr/include/hdf5/serial
  CFLAGS += -DHDF5
  LDLIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
  #LDLIBS += /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a
  #LDLIBS += /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a
  LDLIBS += -lhdf5_hl -lhdf5 -lpthread -ldl -lsz -lz -lm
  #endif
